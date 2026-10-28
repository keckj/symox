"""Differential tests for wave-19 task polycap1.

* ``count_roots`` over ZZ/QQ uses an integer Sturm chain (positive multiples
  of the reference's QQ Sturm sequence, lazily scaled pseudo-remainders and
  homogeneous integer evaluation).  Expected counts are hard-coded from
  SymPy 1.14, including bounds that are roots, reversed bounds, rational,
  Float and infinite bounds.
* The modular factorization charges dense divisions by their actual
  operation count, so sparse trinomials the reference factors are no longer
  refused by the work estimate.  Results are checked against SymPy 1.14
  (same factors, same order) and by exact products.
"""
import pytest

from symox import (Float, Integer, Poly, Rational as R, Symbol, chebyshevt, count_roots,
                   expand, factor, factor_list, oo, quo, I)

x = Symbol('x')

POLYS = {
    'x**400 - x - 1': x**400 - x - 1,
    'x**401 - x - 1': x**401 - x - 1,
    'x**2 - 1': x**2 - 1,
    '(x**2 - 2)**3*(x - 1)': expand((x**2 - 2)**3*(x - 1)),
    'chebyshevt(40, x)': chebyshevt(40, x),
    'R(-3, 7)*x**5 + x/2 - 1': R(-3, 7)*x**5 + x/2 - 1,
    '-x**301 + 5*x**4 - x + R(1, 3)': -x**301 + 5*x**4 - x + R(1, 3),
    'expand((x**3 - x + R(1, 7))*(x**2 - 3)**2)': expand((x**3 - x + R(1, 7))*(x**2 - 3)**2),
    'x**6 - 7*x**4 + 14*x**2 - 7': x**6 - 7*x**4 + 14*x**2 - 7,
}
BOUNDS = [None, (0, 2), (-2, 0), (R(-3, 2), R(7, 5)), (-1.5, Float('1.25')), (1, -1),
          (-1, 1), (R(1, 3), oo), (-oo, R(-1, 2)), (0, 0), (1, 1)]
# SymPy 1.14 count_roots(f, *bounds) for each entry of BOUNDS
EXPECTED = [
    ('x**400 - x - 1', [2, 1, 1, 2, 2, 1, 1, 1, 1, 0, 0]),
    ('x**401 - x - 1', [1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0]),
    ('x**2 - 1', [2, 1, 1, 2, 2, 2, 2, 1, 1, 0, 1]),
    ('(x**2 - 2)**3*(x - 1)', [3, 2, 1, 2, 2, 2, 1, 2, 1, 0, 1]),
    ('chebyshevt(40, x)', [40, 20, 20, 40, 40, 40, 40, 16, 13, 0, 0]),
    ('R(-3, 7)*x**5 + x/2 - 1', [1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0]),
    ('-x**301 + 5*x**4 - x + R(1, 3)', [1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0]),
    ('expand((x**3 - x + R(1, 7))*(x**2 - 3)**2)', [5, 3, 2, 3, 3, 2, 2, 2, 2, 0, 0]),
    ('x**6 - 7*x**4 + 14*x**2 - 7', [6, 3, 3, 2, 2, 2, 2, 3, 3, 0, 0]),
]


@pytest.mark.parametrize('name,counts', EXPECTED)
def test_count_roots_matches_reference(name, counts):
    f = POLYS[name]
    p = Poly(f, x)
    for bounds, expected in zip(BOUNDS, counts):
        args = () if bounds is None else bounds
        got = count_roots(f, *args)
        assert got == expected and isinstance(got, Integer), (bounds, got)
        assert p.count_roots(*args) == expected


def test_count_roots_high_degree_trinomial():
    # The reference's QQ Sturm sequence of x**n - x - 1 ends with a linear
    # element; dividing the derivative by it is the expensive step.
    assert count_roots(x**1200 - x - 1) == 2
    assert count_roots(x**1201 - x - 1) == 1
    assert count_roots(x**1200 - x - 1, 0, 2) == 1
    assert count_roots(x**1200 - x - 1, R(-3, 2), R(-1, 2)) == 1


def test_count_roots_complex_rectangle_unchanged():
    # complex counts keep the rectangle algorithm (SymPy 1.14 values)
    assert count_roots(x**4 + 1, -2 - 2*I, 2 + 2*I) == 4
    assert count_roots(x**4 + 1, 0, 2 + 2*I) == 1
    assert count_roots(x**12 - x - 1, -2 - 2*I, 2 + 2*I) == 12


def _cofactor(f, g):
    q = quo(f, g)
    assert expand(q*g - f) == 0
    return q


@pytest.mark.parametrize('n', [200, 300])
def test_trinomial_factor_list_matches_reference(n):
    f = x**n + x + 1
    got = factor_list(f)
    if n % 3 == 2:
        # SymPy 1.14: (1, [(x**2 + x + 1, 1), (cofactor, 1)])
        assert got == (1, [(x**2 + x + 1, 1), (_cofactor(f, x**2 + x + 1), 1)])
    else:
        # SymPy 1.14: irreducible
        assert got == (1, [(f, 1)])


def test_trinomial_factor_expression():
    assert factor(x**300 + x + 1) == x**300 + x + 1
    f = x**200 + x + 1
    assert factor(f) == (x**2 + x + 1)*_cofactor(f, x**2 + x + 1)


def test_sparse_high_degree_factor_lists():
    # SymPy 1.14: both irreducible over ZZ
    assert factor_list(x**256 - x**5 + 3) == (1, [(x**256 - x**5 + 3, 1)])
    assert factor_list(x**210 + 2*x**3 - 7) == (1, [(x**210 + 2*x**3 - 7, 1)])
