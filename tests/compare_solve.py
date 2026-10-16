"""Differential tests for solve() against SymPy."""

from __future__ import annotations

import sys

import sympy
import symox

CASES = r"""
solve(x-2, x)
solve(2*x+3, x)
solve(x**2-4, x)
solve(x**2+1, x)
solve(x**2-2, x)
solve(x**2+x-1, x)
solve(a*x**2+b*x+c, x)
solve(x**3-6*x**2+11*x-6, x)
solve(x**3-1, x)
solve(x**4-1, x)
solve(Eq(x**2, 9), x)
solve(exp(x)-2, x)
solve(log(x)-1, x)
solve(sin(x), x)
solve(cos(x)-1, x)
solve(sin(x)-Rational(1,2), x)
solve(sqrt(x)-3, x)
solve(1/x-2, x)
solve((x-1)/(x-2), x)
solve(x**2-y, x)
solve(exp(2*x)-3*exp(x)+2, x)
solve([x+y-3, x-y-1], [x, y])
solve([x+y-3, x-y-1], x, y)
solve([x**2+y**2-1, x-y], [x, y])
solve(x*y-1, y)
solve(x**2, x)
solve(x**5-x, x)
solve(2**x-8, x)
solve(tan(x)-1, x)
solve(x**2-2*x+1, x)
solve([2*x+3*y-1, 4*x+6*y-2],[x,y])
solve(x**6-1, x)
solve(x**2-x, x, dict=True)
solve(x**4-5*x**2+4, x)
solve(x**3-2, x)
solve(x - y + 1, y)
solve([x + y + z - 6, x - y - 1, 2*x + z - 5], [x, y, z])
solve(x**2 - 3*x + 1, x)
""".strip().splitlines()


def ns(mod):
    d = {n: getattr(mod, n) for n in dir(mod) if not n.startswith("_")}
    d.update(dict(zip("xyzabc", mod.symbols("x y z a b c"))))
    return d


def run(verbose=True):
    sx, mx = ns(sympy), ns(symox)
    bad, stats = [], {}
    for c in CASES:
        want = eval(c, sx)
        try:
            got = eval(c, mx)
        except Exception as e:  # noqa: BLE001
            got = f"error: {type(e).__name__}: {e}"
        if str(got) == str(want):
            st = "same"
        else:
            # same set of solutions in another order / form?
            try:
                g = sympy.sympify(str(got))
                if isinstance(want, list) and isinstance(g, (list, tuple)) and len(g) == len(want):
                    ok = all(any(sympy.simplify(sympy.sympify(u) - v) == 0 if not isinstance(v, (tuple, dict)) else str(u) == str(v) for u in g) for v in want)
                else:
                    ok = False
            except Exception:  # noqa: BLE001
                ok = False
            st = "equivalent" if ok else "DIFFERENT"
        stats[st] = stats.get(st, 0) + 1
        if st != "same":
            if verbose:
                print(f"[{st:10s}] {c:45s} symox: {str(got):40s} sympy: {want}")
            if st == "DIFFERENT":
                bad.append((c, str(got), str(want)))
    if verbose:
        print("\nsummary:", stats)
    return bad


if __name__ == "__main__":
    sys.exit(1 if run() else 0)
