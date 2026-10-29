"""Dense ZZ subresultant PRS behind resultant(..., includePRS=True).

The sequences are hard-coded from SymPy 1.14 (dup_inner_subresultants).
The resultant values equal SymPy 1.14 except the documented sign correction
(first operand of smaller degree, both degrees odd: the Sylvester
determinant is returned, as before this kernel).
"""
import hashlib
import sys

from symox import (Poly, Rational, expand, resultant, subresultants, symbols,
                   sympify)

x = symbols('x')
R = Rational

CASES = [
    ('x**4 + 1', 'x**2 + 1', '4',
     ['x**4 + 1', 'x**2 + 1', '-2']),
    ('x - 1', 'x**3 + 2', '3',
     ['x**3 + 2', 'x - 1', '-3']),
    ('x**3 + 2', 'x - 1', '-3',
     ['x**3 + 2', 'x - 1', '-3']),
    ('x**3 - 2*x', 'x**5 + 3', '-69',
     ['x**5 + 3', 'x**3 - 2*x', '-4*x - 3', '69']),
    ('2*x**3 + 3*x - 7', '4*x**2 - 6', '-320',
     ['2*x**3 + 3*x - 7', '4*x**2 - 6', '96*x - 112', '-320']),
    ('x**8 + x**6 - 3*x**4 - 3*x**3 + 8*x**2 + 2*x - 5', '3*x**6 + 5*x**4 - 4*x**2 - 9*x + 21', '260708',
     ['x**8 + x**6 - 3*x**4 - 3*x**3 + 8*x**2 + 2*x - 5', '3*x**6 + 5*x**4 - 4*x**2 - 9*x + 21', '15*x**4 - 3*x**2 + 9', '65*x**2 + 125*x - 245', '9326*x - 12300', '260708']),
    ('x**6 - 1', 'x**3 - 1', '0',
     ['x**6 - 1', 'x**3 - 1']),
    ('(x - 1)**3*(x + 2)', '(x - 1)**2*(x**2 + 5)', '0',
     ['x**4 - x**3 - 3*x**2 + 5*x - 2', 'x**4 - 2*x**3 + 6*x**2 - 10*x + 5', '-x**3 + 9*x**2 - 15*x + 7', '54*x**2 - 108*x + 54']),
    ('7*x**5 - 3', '2*x**2', '288',
     ['7*x**5 - 3', '2*x**2', '-48']),
    ('x**10 - 2*x**3 + 1', 'x**7 + 5*x - 11', '7292213145',
     ['x**10 - 2*x**3 + 1', 'x**7 + 5*x - 11', '-5*x**4 + 9*x**3 + 1', '6686*x**3 + 225*x**2 + 3530*x - 6146', '166745*x**2 - 674954*x + 674314', '54081876*x - 72460821', '7292213145']),
]


def test_small_prs_match_the_reference():
    for a, b, value, prs in CASES:
        f, g = expand(sympify(a)), expand(sympify(b))
        r, sequence = resultant(f, g, x, includePRS=True)
        assert str(r) == value, (a, b)
        assert [str(p) for p in sequence] == prs, (a, b)
        assert [str(p) for p in subresultants(f, g, x)] == prs, (a, b)


def test_poly_method_and_other_domains():
    value, prs = Poly(x**5 + 3, x).resultant(Poly(2*x**2 - x, x), includePRS=True)
    assert value == 291
    assert prs == [Poly(x**5 + 3, x), Poly(2*x**2 - x, x), Poly(x + 48, x), Poly(291, x)]
    assert subresultants(x**5 + 3, 2*x**2 - x, x) == [x**5 + 3, 2*x**2 - x, x + 48, 291]
    # QQ keeps the sparse PRS
    assert resultant(x/2 + 1, x**2 + 3, x, includePRS=True) == (
        R(7, 4), [x**2 + 3, x/2 + 1, R(7, 4)])
    value, prs = resultant(Poly(x**3 + 2*x + 1, x, domain='QQ'),
                           Poly(x**2 - 3, x, domain='QQ'), includePRS=True)
    assert value == -74
    assert prs == [Poly(x**3 + 2*x + 1, x, domain='QQ'), Poly(x**2 - 3, x, domain='QQ'),
                   Poly(5*x + 1, x, domain='QQ'), Poly(-74, x, domain='QQ')]
    assert all(p.domain == Poly(1, x, domain='QQ').domain for p in prs)


def _digest(obj):
    old = sys.get_int_max_str_digits()
    sys.set_int_max_str_digits(0)
    try:
        return hashlib.sha256(str(obj).encode()).hexdigest()
    finally:
        sys.set_int_max_str_digits(old)


def test_growing_prs_match_the_reference():
    expected = {
        (40, 30): (32, '49492a4722deb1421aa18674e5b2c924537469e823af6c17b80a96a6b190aa32',
                   '7ea43cfd87d331e670e16cd7d709fd35fa1e61f1bdb027ff9a8a6d23b95347df'),
        (60, 45): (47, '0b629bbc789f1f8305a54270d5020e351df1f667983bce68e85dd5d5d4b319df',
                   'cf5f02f9c418d8206b89430fe192b14c2595e83720ba2c994ef371486c01ed7c'),
    }
    for (n, m), (length, prs_digest, value_digest) in expected.items():
        f, g = expand((x + 1)**n), expand((x - 2)**m + 1)
        r, prs = resultant(f, g, x, includePRS=True)
        assert len(prs) == length
        assert _digest(prs) == prs_digest
        assert _digest(r) == value_digest
        assert r == ((-3)**m + 1)**n
        assert _digest(subresultants(f, g, x)) == prs_digest
