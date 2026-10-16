"""Differential tests for assumptions (Symbol('x', positive=True), ...)."""

from __future__ import annotations

import sys

import sympy
import symox

CASES = r"""
sqrt(p**2)
sqrt(r**2)
sqrt(x**2)
(p**2)**Rational(1,3)
(p**3)**y
sqrt(n**2)
Abs(p)
Abs(n)
Abs(r)
Abs(r)**2
Abs(x)**2
Abs(nn)
sign(p)
sign(n)
sqrt(p*q)
sqrt(p*x)
(p*x)**y
(n*x)**Rational(1,2)
log(exp(r))
log(exp(x))
exp(log(p))
floor(k)
ceiling(k)
sin(k*pi)
cos(k*pi)
tan(k*pi)
(-1)**(2*k)
(-1)**(2*k+1)
re(r)
im(r)
conjugate(r)
p > 0
p < 0
r**2 >= 0
n + p > 0
p + q > 0
n - p < 0
Max(p, n)
Min(p, 0)
Max(p, 0)
sqrt(-p)
(-p)**Rational(1,3)
Abs(p*x)
Abs(I*r)
sqrt(p)**2
(p**Rational(1,2))**Rational(1,2)
(x**2)**y
(-p)**k
p.is_positive
n.is_positive
r.is_positive
(p + 1).is_positive
(p*n).is_negative
(r**2).is_nonnegative
(r**2 + 1).is_positive
exp(r).is_positive
x.is_real
(r + p).is_real
(2*k + 1).is_integer
(2*k + 1).is_odd
(2*k).is_even
sqrt(p).is_positive
(pi - 3).is_positive
I.is_real
k.is_integer
(k/2).is_integer
srepr(p)
srepr(r)
p == Symbol('p', positive=True)
p == Symbol('p')
""".strip().splitlines()


def ns(mod):
    d = {n: getattr(mod, n) for n in dir(mod) if not n.startswith("_")}
    S = mod.Symbol
    d.update(
        x=S("x"), y=S("y"),
        p=S("p", positive=True), q=S("q", positive=True), n=S("n", negative=True),
        r=S("r", real=True), k=S("k", integer=True), nn=S("nn", nonnegative=True),
    )
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
                print(f"[{st}] {c:32s} symox: {str(got):30s} sympy: {want}")
            bad.append((c, str(got), str(want)))
    if verbose:
        print("\nsummary:", stats)
    return bad


if __name__ == "__main__":
    sys.exit(1 if run() else 0)
