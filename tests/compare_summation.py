"""Differential tests for summation() against SymPy."""

from __future__ import annotations

import sys

import sympy
import symox

CASES = r"""
summation(k, (k, 1, n))
summation(k**2, (k, 1, n))
summation(k**3, (k, 0, n))
summation(k**4, (k, 1, n))
summation(2*k + 1, (k, 0, n - 1))
summation(k**2, (k, m, n))
summation(x, (k, 1, n))
summation(k, (k, 1, 10))
summation(k**2 + 3*k, (k, 1, 100))
summation(1/k, (k, 1, 5))
summation(2**k, (k, 0, n))
summation(3**k, (k, 1, n))
summation(x**k, (k, 0, n))
summation(k*2**k, (k, 0, n))
summation(Rational(1, 2)**k, (k, 0, oo))
summation(Rational(1, 3)**k, (k, 1, oo))
summation(k*Rational(1, 2)**k, (k, 0, oo))
summation(1/k**2, (k, 1, oo))
summation(1/k**4, (k, 1, oo))
summation(1/k**6, (k, 1, oo))
summation(x*k**2, (k, 1, n))
summation(k*(k + 1), (k, 1, n))
""".strip().splitlines()


def ns(mod):
    d = {n: getattr(mod, n) for n in dir(mod) if not n.startswith("_")}
    d.update(dict(zip(["k", "n", "m", "x"], mod.symbols("k n m x"))))
    return d


def run(verbose=True):
    sx, mx = ns(sympy), ns(symox)
    bad, stats = [], {}
    for c in CASES:
        want = eval(c, sx)
        try:
            got = eval(c, mx)
        except Exception as e:  # noqa: BLE001
            got = f"error: {e}"
        if str(got) == str(want):
            st = "same"
        else:
            try:
                back = sympy.sympify(str(got))
                diff = sympy.simplify(back - want)
                st = "equivalent" if diff == 0 else "DIFFERENT"
            except Exception:  # noqa: BLE001
                st = "DIFFERENT"
        stats[st] = stats.get(st, 0) + 1
        if st != "same" and verbose:
            print(f"[{st:10s}] {c:40s} symox: {str(got):40s} sympy: {want}")
        if st == "DIFFERENT":
            bad.append((c, str(got), str(want)))
    if verbose:
        print("\nsummary:", stats)
    return bad


if __name__ == "__main__":
    sys.exit(1 if run() else 0)
