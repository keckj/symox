"""Differential tests for integrate() against SymPy.

An antiderivative is accepted when it prints exactly like SymPy's, or when
its derivative equals the integrand numerically (antiderivatives may differ
by a constant). Definite integrals must agree numerically.
"""

from __future__ import annotations

import cmath
import random
import sys

import sympy
import symox

INDEFINITE = r"""
x
x**2
x**5
1/x
1/x**2
sqrt(x)
x**Rational(3,2)
exp(x)
exp(2*x)
exp(-3*x+1)
sin(x)
cos(x)
tan(x)
sin(2*x+1)
cos(3*x)
x*exp(x)
x**2*exp(x)
x*sin(x)
x*cos(x)
x**2*sin(x)
log(x)
x*log(x)
x**3*log(x)
log(x)**2
atan(x)
asin(x)
x*atan(x)
1/(1+x**2)
1/(4+x**2)
1/sqrt(1-x**2)
1/sqrt(1+x**2)
sqrt(1-x**2)
sin(x)**2
cos(x)**2
sin(x)**3
cos(x)**3
sin(x)*cos(x)
sin(x)**2*cos(x)**2
sin(x)**3*cos(x)**2
sec(x)**2
tan(x)**2
exp(x)*sin(x)
exp(x)*cos(x)
exp(2*x)*sin(3*x)
2*x*exp(x**2)
x/(x**2+1)
cos(x)*exp(sin(x))
sinh(x)
cosh(x)
tanh(x)
1/(x+1)
1/(2*x+3)
1/(2*x+3)**3
x/(x+1)**2
(x**2+1)/(x-1)
1/(x**2-1)
1/(x**2+2*x+5)
(x+1)/(x**2+x+1)
(3*x+2)/(x**2-5*x+6)
3*x**2+2*x+1
(2*x+3)**5
exp(x)/(1+exp(x))
log(x)/x
1/(x*log(x))
exp(-x**2)
2**x
x*y
y*exp(y*x)
x**y
sin(y*x)
""".strip().splitlines()

DEFINITE = r"""
x**2, (x, 0, 1)
sin(x), (x, 0, pi)
exp(-x), (x, 0, 1)
1/(1+x**2), (x, 0, 1)
x*exp(x), (x, 0, 2)
log(x), (x, 1, E)
cos(x)**2, (x, 0, pi)
""".strip().splitlines()


def ns(mod):
    d = {n: getattr(mod, n) for n in dir(mod) if not n.startswith("_")}
    d.update(dict(zip("xyz", mod.symbols("x y z"))))
    return d


def num_check(expr_s, f_s, syms, deriv=True):
    rnd = random.Random(7)
    d = sympy.diff(expr_s, syms[0]) if deriv else expr_s
    fd = sympy.lambdify(syms, d, "mpmath")
    ff = sympy.lambdify(syms, f_s, "mpmath")
    ok = 0
    for _ in range(6):
        pt = [complex(rnd.uniform(0.3, 0.9), rnd.uniform(-0.2, 0.2)) for _ in syms]
        try:
            a, b = complex(fd(*pt)), complex(ff(*pt))
        except Exception:  # noqa: BLE001
            continue
        if abs(a - b) > 1e-7 * max(1, abs(b)):
            return False
        ok += 1
    return ok > 0


def run(verbose=True):
    bad = []
    stats = {}
    sx = ns(sympy)
    mx = ns(symox)
    for case in INDEFINITE:
        want = sympy.integrate(eval(case, sx), sx["x"])
        try:
            got = symox.integrate(eval(case, mx), mx["x"])
        except Exception as e:  # noqa: BLE001
            got = f"error {e}"
        if str(got) == str(want):
            status = "same"
        elif "Integral" in str(got):
            status = "unevaluated"
        else:
            f_s = eval(case, sx)
            back = sympy.sympify(str(got))
            syms = [sx["x"]] + sorted((f_s.free_symbols | back.free_symbols) - {sx["x"]}, key=str)
            status = "verified" if num_check(back, f_s, syms) else "WRONG"
        stats[status] = stats.get(status, 0) + 1
        if status != "same" and verbose:
            print(f"[{status:11s}] {case:28s} symox: {str(got):40s} sympy: {want}")
        if status in ("WRONG", "unevaluated"):
            bad.append((case, str(got), str(want)))
    for case in DEFINITE:
        f, lim = case.split(", ", 1)
        want = sympy.integrate(eval(f, sx), eval(lim, sx))
        got = symox.integrate(eval(f, mx), eval(lim, mx))
        same = str(got) == str(want)
        ok = same or abs(complex(sympy.N(sympy.sympify(str(got)))) - complex(sympy.N(want))) < 1e-10
        status = "same" if same else ("verified" if ok else "WRONG")
        stats[status] = stats.get(status, 0) + 1
        if status != "same" and verbose:
            print(f"[{status:11s}] {case:28s} symox: {str(got):40s} sympy: {want}")
        if not ok:
            bad.append((case, str(got), str(want)))
    if verbose:
        print("\nsummary:", stats)
    return bad


if __name__ == "__main__":
    sys.exit(1 if run() else 0)
