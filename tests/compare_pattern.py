"""Differential tests for pattern matching (Wild, match, replace)."""

from __future__ import annotations

import sys

import sympy
import symox

CASES = r"""
(2*x+3).match(a*x+b)
x.match(a*x+b)
(2*x+3).match(ax*x+bx)
(x**2).match(x**a)
x.match(x**a)
(sin(x)**2).match(sin(a)**b)
(3*x*y).match(a*x)
(3*x*y).match(ax*x)
(exp(2*x)).match(exp(ax*x))
(2*x**3+x).match(ax*x**3+bx*x)
(y*x+z).match(ax*x+bx)
(x+y).match(x+a)
x.match(y)
(2*x).match(a)
y.match(ax*x)
(sin(x)+cos(x)).replace(sin(a), cos(a))
(sin(x)**2+sin(y)).replace(sin(a), a)
(f(x)+f(y)).replace(f(a), a**2)
(sin(x)+sin(y)).replace(sin, cos)
(sin(x)*cos(x)).replace(sin, lambda t: t**2)
(x**2 + y**2).replace(x**2, z)
(sin(x)**2 + cos(x)).replace(sin(a)**2, lambda a: 1 - cos(a)**2)
""".strip().splitlines()


def ns(mod):
    d = {n: getattr(mod, n) for n in dir(mod) if not n.startswith("_")}
    x, y, z = mod.symbols("x y z")
    d.update(x=x, y=y, z=z, f=mod.Function("f"))
    d.update(a=mod.Wild("a"), b=mod.Wild("b"), ax=mod.Wild("ax", exclude=[x]), bx=mod.Wild("bx", exclude=[x]))
    return d


def canon(v):
    if isinstance(v, dict):
        return sorted((str(k), str(w)) for k, w in v.items())
    return str(v)


def run(verbose=True):
    sx, mx = ns(sympy), ns(symox)
    bad, stats = [], {}
    for c in CASES:
        want = eval(c, sx)
        try:
            got = eval(c, mx)
        except Exception as e:  # noqa: BLE001
            got = f"error: {type(e).__name__}: {e}"
        st = "same" if canon(got) == canon(want) else "DIFFERENT"
        stats[st] = stats.get(st, 0) + 1
        if st != "same":
            if verbose:
                print(f"[{st}] {c:50s} symox: {canon(got)} sympy: {canon(want)}")
            bad.append((c, str(got), str(want)))
    if verbose:
        print("\nsummary:", stats)
    return bad


if __name__ == "__main__":
    sys.exit(1 if run() else 0)
