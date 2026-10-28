"""Exact evaluation of the large combinatorial values (TODO section 8).

The expected digests are sha256(str(value))[:16] of the SymPy 1.14 results
(pure Python ground types), recorded by the harness in
audits/performance/2026-10-29-bigexact/.
"""
import hashlib

import pytest

from symox import (
    Integer, Mod, Rational, bernoulli, binomial, catalan, euler, factorial, factorial2,
    ff, fibonacci, gamma, genocchi, legendre, lucas, primorial, rf,
    subfactorial, tribonacci, zoo,
)


def digest(value):
    return hashlib.sha256(text_of(value).encode()).hexdigest()[:16]


def text_of(value):
    # a Python int above 4300 digits has no str() by default
    return str(Integer(value)) if isinstance(value, int) else str(value)


SYMPY_DIGESTS = [
    # the TODO witnesses, previously left unevaluated
    (lambda: factorial(100001), 456579, 'e5ac91c2bf24cae9'),
    (lambda: fibonacci(1000001), 208988, '587ad4a35a5b4069'),
    (lambda: catalan(100001), 60199, '878ea0c8848d7e3a'),
    (lambda: bernoulli(5002), 12356, '73638452cad54074'),
    (lambda: legendre(2001, Rational(1, 3)), 1907, '38fd57a29de72222'),
    # neighbours sharing the caps
    (lambda: gamma(100002), 456579, 'e5ac91c2bf24cae9'),
    (lambda: lucas(1000001), 208988, '8c5415c2f2500fbc'),
    (lambda: tribonacci(100001), 26465, '3b83e22a2708f24d'),
    (lambda: rf(3, 100001), 456589, 'e09e51f019c4a699'),
    (lambda: factorial2(200001), 486680, 'b73bf179b26aad40'),
    (lambda: euler(5002), 15354, '721b78e5d18fea59'),
    (lambda: genocchi(5002), 13848, '85dc43613a07998a'),
    # smaller instances of the same code paths
    (lambda: factorial(3001), 9135, '4b4c88fbd476f646'),
    (lambda: fibonacci(20001), 4180, '2465e7bd7ea89ee8'),
    (lambda: catalan(2001), 1200, '1a774cdf20c4a711'),
    (lambda: bernoulli(700), 1155, 'a048d54155b265a4'),
    (lambda: bernoulli(1202), 2226, '01224a0bb9e741dd'),
    (lambda: legendre(301, Rational(-2, 7)), 685, '5dae99be813ecea9'),
    (lambda: legendre(150, Rational(5, 4)), 315, '22f56ae266d9e27a'),
    (lambda: lucas(20001), 4180, 'af97c110abca809a'),
    (lambda: subfactorial(1001), 2571, '3e0c7feaffa50f79'),
    (lambda: tribonacci(3001), 794, 'bf0d0548668a7695'),
    (lambda: ff(10**6, 1001), 6006, '5dfbaad19402060e'),
    (lambda: rf(-3000, 1001), 3400, 'a49397f53cd88f5b'),
    (lambda: ff(-7, 300), 627, 'a140de0e779a4ec2'),
    (lambda: factorial2(3001), 4568, '075c32ec755f4a33'),
    (lambda: factorial2(3000), 4567, '403e3b99440823a7'),
    (lambda: euler(700), 1553, '73afceea84d69141'),
    (lambda: euler(1202), 2948, '27d8aa710bccfe13'),
    (lambda: genocchi(702), 1348, '9310866421ee6763'),
    (lambda: binomial(3001, 1000), 828, 'f07bdc30acd4802a'),
    (lambda: binomial(-3001, 1000), 976, '4420f616b00ad487'),
    (lambda: primorial(3000), 11829, '41f2f0b447083028'),
]


@pytest.mark.parametrize('case', range(len(SYMPY_DIGESTS)))
def test_matches_sympy_digest(case):
    make, length, expected = SYMPY_DIGESTS[case]
    value = make()
    # primorial returns a Python int, as in SymPy
    assert isinstance(value, int) or value.is_Rational
    text = text_of(value)
    assert (len(text), digest(value)) == (length, expected)


def test_small_exact_values():
    assert legendre(-41, Rational(1, 3)) == Rational(571570263174244609, 12157665459056928801)
    assert legendre(41, Rational(1, 3)) == legendre(-42, Rational(1, 3))
    assert legendre(7, 0) == 0 and legendre(8, 0) == Rational(35, 128)
    assert legendre(9, -1) == -1 and legendre(6, 2) == Rational(10159, 16)
    assert rf(5, -300) is zoo
    assert rf(-5, 3) == -60 and ff(-5, 3) == -210 and ff(5, -2) == Rational(1, 42)
    assert subfactorial(0) == 1 and subfactorial(1) == 0 and subfactorial(10) == 1334961
    assert tribonacci(0) == 0 and tribonacci(3) == 2 and tribonacci(10) == 149
    assert euler(200001) == 0 and bernoulli(1000001) == 0


def test_catalan_and_legendre_certificates():
    n = 3000
    c = catalan(n)
    assert c*(n + 1) == binomial(2*n, n) == factorial(2*n)/factorial(n)**2
    # P_n(x) through the explicit sum 2^-n sum C(n,k)^2 (x-1)^(n-k) (x+1)^k
    x = Rational(3, 7)
    n = 60
    explicit = sum(binomial(n, k)**2*(x - 1)**(n - k)*(x + 1)**k for k in range(n + 1))/2**n
    assert legendre(n, x) == explicit


def test_enormous_values_stay_unevaluated():
    for expr in [factorial(10**9), fibonacci(10**12), lucas(10**12), catalan(10**9),
                 tribonacci(10**12), subfactorial(10**9), bernoulli(10**7),
                 euler(10**7), binomial(10**12, 10**9), rf(3, 10**9),
                 legendre(10**7, Rational(1, 3))]:
        assert not expr.is_Number
    # an unevaluated huge integer never takes an inexact Mod
    m = Mod(fibonacci(10**12), 1000)
    assert m.func is Mod or m.is_Integer
