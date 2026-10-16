"""Differential tests for the targeted simplification routines."""

from __future__ import annotations

import sys

import sympy
import symox

CASES = r"""
expand_trig(sin(2*x))
expand_trig(cos(2*x))
expand_trig(sin(x+y))
expand_trig(cos(x+y))
expand_trig(sin(3*x))
expand_trig(cos(3*x))
trigsimp(sin(x)**2+cos(x)**2)
trigsimp(2*sin(x)**2+2*cos(x)**2)
trigsimp(sin(x)/cos(x))
trigsimp(1-sin(x)**2)
trigsimp(2*sin(x)*cos(x))
trigsimp(cos(x)**2-sin(x)**2)
trigsimp(sin(x)**4 - cos(x)**4)
trigsimp(1+tan(x)**2)
trigsimp(sin(x)*cos(y)+cos(x)*sin(y))
trigsimp(cosh(x)**2-sinh(x)**2)
trigsimp(x*sin(y)**2 + x*cos(y)**2)
powsimp(x**a*x**b)
powsimp(x**a*y**a)
powsimp(exp(x)*exp(y))
powsimp(2**x*2**y)
powsimp(x**2*x**a)
expand_log(log(p*q))
expand_log(log(p**2))
expand_log(log(x*y))
expand_log(log(x*y), force=True)
logcombine(log(p)+log(q))
logcombine(2*log(p))
logcombine(log(x)+log(y))
logcombine(log(x)+log(y), force=True)
expand_power_base((x*y)**a, force=True)
expand_power_exp(x**(a+b))
expand_mul(x*(y+z))
expand(sin(2*x), trig=True)
expand(log(p*q))
simplify(sin(x)**2+cos(x)**2)
simplify(2*sin(x)*cos(x))
nsimplify(0.5)
nsimplify(0.333333333333333)
nsimplify(1.4142135623731)
nsimplify(3.14159265358979, [pi])
nsimplify(0.25*x + 1.5)
""".strip().splitlines()


def ns(mod):
    d = {n: getattr(mod, n) for n in dir(mod) if not n.startswith("_")}
    x, y, z, a, b = mod.symbols("x y z a b")
    p, q = mod.symbols("p q", positive=True)
    d.update(x=x, y=y, z=z, a=a, b=b, p=p, q=q)
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
        st = "same" if str(got) == str(want) else "DIFFERENT"
        stats[st] = stats.get(st, 0) + 1
        if st != "same":
            if verbose:
                print(f"[{st}] {c:45s} symox: {str(got):30s} sympy: {want}")
            bad.append((c, str(got), str(want)))
    if verbose:
        print("\nsummary:", stats)
    return bad


if __name__ == "__main__":
    sys.exit(1 if run() else 0)
