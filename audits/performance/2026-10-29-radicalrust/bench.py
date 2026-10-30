"""Timing harness for the native unrad, radical and negative-base steps.

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
    "unrad(sqrt(x) + sqrt(x + 1) + sqrt(1 - x))",
    "unrad(root(x, 3) - root(x + 1, 4)/2 + root(x + 2, 3))",
    "unrad(sqrt(x) + root(x, 3) + y)",
    "solve(sqrt(x) + sqrt(x + 1) - 3, x)",
    "solve(sqrt(x) + root(x, 3) - 2, x)",
    "solve(root(x + 1, 3) + root(x, 3) - 1, x)",
    "solveset(sqrt(2*x - 1) - sqrt(x - 4) - 2, x, S.Reals)",
    "solveset(-8*x**2/(9*(x**2 - 1)**(S(4)/3)) + 4/(3*(x**2 - 1)**(S(1)/3)), x, S.Complexes)",
    "solveset((-2)**x + 8, x, S.Reals)",
    "solveset(3*(-2)**x - 12, x, S.Reals)",
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
    from importlib import import_module
    ns = {k: getattr(M, k) for k in ('solve', 'solveset', 'sqrt', 'root', 'S')}
    ns['unrad'] = import_module(M.__name__ + '.solvers.solvers').unrad
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
