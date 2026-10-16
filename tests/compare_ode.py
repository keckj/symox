"""Differential tests for dsolve() against SymPy.

A solution is accepted when it prints like SymPy's or when it satisfies
the ODE (checked by SymPy's checkodesol) with the same number of
arbitrary constants.
"""

from __future__ import annotations

import sys

import sympy
import symox

CASES = r"""
dsolve(f(x).diff(x) - f(x), f(x))
dsolve(f(x).diff(x) + 2*f(x), f(x))
dsolve(f(x).diff(x) - x, f(x))
dsolve(f(x).diff(x) - x*f(x), f(x))
dsolve(f(x).diff(x) + f(x) - x, f(x))
dsolve(f(x).diff(x, 2) - f(x), f(x))
dsolve(f(x).diff(x, 2) + f(x), f(x))
dsolve(f(x).diff(x, 2) + 2*f(x).diff(x) + f(x), f(x))
dsolve(f(x).diff(x, 2) - 3*f(x).diff(x) + 2*f(x), f(x))
dsolve(f(x).diff(x, 2) + 4*f(x) - sin(x), f(x))
dsolve(f(x).diff(x, 2) - f(x) - exp(2*x), f(x))
dsolve(f(x).diff(x) - f(x)**2, f(x))
dsolve(f(x).diff(x) - a*f(x), f(x))
dsolve(f(x).diff(x,3) - f(x), f(x))
dsolve(Eq(f(x).diff(x), f(x)), f(x), ics={f(0): 1})
dsolve(f(x).diff(x, 2) + f(x), f(x), ics={f(0): 0, f(x).diff(x).subs(x, 0): 1})
dsolve(f(x).diff(x) - f(x)/x, f(x))
dsolve(x**2*f(x).diff(x,2) + x*f(x).diff(x) - f(x), f(x))
dsolve(f(x).diff(x, 2) - 2*f(x).diff(x) + 5*f(x), f(x))
dsolve(f(x).diff(x) - sin(x)*f(x), f(x))
dsolve(f(x).diff(x, 2) - x, f(x))
dsolve(f(x).diff(x, 2) - f(x) - exp(x), f(x))
dsolve(f(x).diff(x, 4) - f(x), f(x))
dsolve(f(x).diff(x, 2) + 9*f(x) - cos(3*x), f(x))
""".strip().splitlines()


def ns(mod):
    d = {n: getattr(mod, n) for n in dir(mod) if not n.startswith("_")}
    d.update(x=mod.Symbol("x"), a=mod.Symbol("a"), f=mod.Function("f"))
    return d


def nconst(e):
    return len({s for s in e.free_symbols if str(s).startswith("C")})


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
            try:
                g = sympy.sympify(str(got), locals={"f": sympy.Function("f")})
                ode = eval(c.split(", f(x)")[0].replace("dsolve(", "", 1), sx)
                ok = sympy.checkodesol(ode, g)[0] and nconst(g) == nconst(want)
            except Exception as e:  # noqa: BLE001
                ok = False
            st = "verified" if ok else "DIFFERENT"
        stats[st] = stats.get(st, 0) + 1
        if st != "same":
            if verbose:
                print(f"[{st:9s}] {c:60s}\n            symox: {got}\n            sympy: {want}")
            if st == "DIFFERENT":
                bad.append((c, str(got), str(want)))
    if verbose:
        print("\nsummary:", stats)
    return bad


if __name__ == "__main__":
    sys.exit(1 if run() else 0)
