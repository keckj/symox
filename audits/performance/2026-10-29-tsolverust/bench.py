"""Timing harness for the native _tsolve steps.

Usage: python bench.py (symox|sympy) (first|repeat) [repeats]
  first : one call per case, each case in a fresh process is driven by run.sh
          (here: times the first call of every case in this process order
          given by --case)
  repeat: median of `repeats` calls of each case in one process (after the
          first call, caches warm).
Prints one JSON line per case.
"""
import json
import statistics
import sys
import time

CASES = [
    "solve(x**y - 2, x)",
    "solve(-a*x + 2*x*log(x), x)",
    "solve(x**(2*y) - 3, x)",
    "solve((x**2 - 1)**y - 4, x)",
    "solve(x*exp(x**2) - 1, x)",
    "solve(sqrt(x)*exp(sqrt(x)) - 1, x)",
    "solve(x - exp(-x**2), x)",
    "solve(x**a*exp(x) - 1, x)",
    "solve(x**a - exp(x), x)",
    "solve(sqrt(x) + log(x), x)",
]


def main():
    lib, mode = sys.argv[1], sys.argv[2]
    reps = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    only = int(sys.argv[4]) if len(sys.argv) > 4 else None
    t0 = time.perf_counter()
    if lib == 'sympy':
        import sympy as M
    else:
        import symox as M
    imp = time.perf_counter() - t0
    ns = {k: getattr(M, k) for k in ('solve', 'log', 'exp', 'sqrt')}
    ns.update(zip('xya', M.symbols('x y a')))
    for i, c in enumerate(CASES):
        if only is not None and i != only:
            continue
        code = compile(c, c, 'eval')
        t = time.perf_counter()
        r = eval(code, ns)
        first = time.perf_counter() - t
        rec = {'case': i, 'expr': c, 'first_s': first, 'import_s': imp,
               'result': str(r)}
        if mode == 'repeat':
            ts = []
            for _ in range(reps):
                t = time.perf_counter()
                eval(code, ns)
                ts.append(time.perf_counter() - t)
            rec['repeat_median_s'] = statistics.median(ts)
            rec['repeats'] = reps
        print(json.dumps(rec), flush=True)


main()
