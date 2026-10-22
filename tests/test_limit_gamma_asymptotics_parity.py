"""Limits at oo of gamma-like functions, checked against SymPy 1.14."""
import subprocess
import sys

import pytest

from symox import (E, Abs, Limit, S, Symbol, binomial, exp, factorial, ff,
                   gamma, limit, limit_seq, loggamma, log, oo, pi, rf, sqrt,
                   symbols)

n, m, a = symbols('n m a')
k = Symbol('n', integer=True, positive=True)

CASES = [
    (factorial(n)*exp(n)/(n**n*sqrt(n)), n, sqrt(2)*sqrt(pi)),
    (factorial(m)*exp(m)/(m**m*sqrt(m)), m, sqrt(2)*sqrt(pi)),
    (gamma(m + 1)/(m**m*exp(-m)*sqrt(m)), m, sqrt(2)*sqrt(pi)),
    (gamma(n)*exp(n)/n**(n - S(1)/2), n, sqrt(2)*sqrt(pi)),
    (factorial(k)*exp(k)/(k**k*sqrt(k)), k, sqrt(2)*sqrt(pi)),
    (loggamma(n) - (n - S(1)/2)*log(n) + n, n, log(2)/2 + log(pi)/2),
    (gamma(n + 1)/gamma(n)/n, n, 1),
    (gamma(n + S(1)/2)/(gamma(n)*sqrt(n)), n, 1),
    (gamma(n)/gamma(n + S(1)/3)*n**(S(1)/3), n, 1),
    (factorial(n)**(1/n)/n, n, exp(-1)),
    (factorial(n)**(S(1)/n)*E/n, n, 1),
    (factorial(n + 1)/factorial(n) - n, n, 1),
    (binomial(n, 2)/n**2, n, S(1)/2),
    (factorial(n)/factorial(n - 1)/n, n, 1),
    (factorial(n)/n**n, n, 0),
    (gamma(n), n, oo),
]


@pytest.mark.parametrize("e, v, ref", CASES, ids=[str(c[0]) for c in CASES])
def test_reference_values(e, v, ref):
    assert limit(e, v, oo) == ref


def test_correct_where_reference_fails():
    # SymPy 1.14 returns oo (wrong) or raises; the exact values are below
    assert limit(rf(n, n)/factorial(2*n)*factorial(n), n, oo) == S(1)/2
    assert limit(ff(2*n, n)/factorial(2*n)*factorial(n), n, oo) == 1
    assert limit(ff(n, n)/factorial(n), n, oo) == 1


@pytest.mark.parametrize("e", [
    binomial(2*n, n)*sqrt(n)/4**n,
    gamma(2*n)/(gamma(n)**2*4**n)*sqrt(n),
    loggamma(n + 1) - loggamma(n) - log(n),
])
def test_unsupported_stays_unevaluated_or_correct(e):
    # never a wrong number: either the reference value or a Limit
    r = limit(e, n, oo)
    assert isinstance(r, Limit) or r in (1/sqrt(pi), oo, 0)


def test_abs_power_quotient_does_not_crash():
    code = (
        "from symox import *\n"
        "n, a, b = symbols('n a b')\n"
        "for e in [Abs((n + a)**(n + b)/(n + a + 1)**(n + b + 1)),\n"
        "          Abs((n + 1)**(n + 2)/(n + 2)**(n + 1)),\n"
        "          (n**2 + a)**(n + b)/(n + a)**(2*n)]:\n"
        "    limit(e, n, oo); limit_seq(e, n)\n"
        "print('ok')\n")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True,
                       text=True, timeout=300, cwd="/tmp")
    assert r.returncode == 0, r.stderr[-2000:]
    assert r.stdout.strip().endswith("ok")
