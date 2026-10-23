"""Solve leftovers of wave 6 in the form of SymPy 1.14.

Expected values are hard-coded from SymPy 1.14:

* ``_tsolve``'s bivariate block (``_solve_lambert``, ``bivariate_type``),
  ported in ``python/symox/solvers/bivariate.py``;
* solveset radical products (``unrad`` on ``sqrt(x)*sqrt(x - 7) - 12``);
* the ``Complement(..., poles)`` of a rational equation with symbolic poles;
* realness of Cardano roots decided after ``expand_complex``;
* ``roots`` with ``preprocess_roots``'s ``_integer_basis`` rescaling;
* limits of differences of log-singular functions are never ``nan``.
"""

import pytest

from symox import (Ei, Limit, N, Rational as R, S, Symbol, exp, limit,
                   log, roots, solve, solveset, sqrt, symbols)

x, y, a, b = symbols('x y a b')


LAMBERT = [
    (lambda: solve(x*exp(x**2) - 1, x), ['exp(-LambertW(2)/2)']),
    (lambda: solve(x**a*exp(x) - 1, x), ['a*LambertW(1/a)']),
    (lambda: solve(log(x) + x**2, x), ['exp(-LambertW(2)/2)']),
    (lambda: solve(x**x - 2, x), ['exp(LambertW(log(2)))']),
    (lambda: solve(x*log(x) - 1, x), ['exp(LambertW(1))']),
    (lambda: solve(x**2 - 2**x, x),
     ['2', '4', '-2*LambertW(log(2)/2)/log(2)']),
    (lambda: solve(exp(x) + x**3, x),
     ['-3*LambertW(1/3)', '-3*LambertW(-1/6 - sqrt(3)*I/6)',
      '-3*LambertW(-1/6 + sqrt(3)*I/6)']),
    (lambda: solve(x - log(x)*y, x), ['-y*LambertW(-1/y)']),
]


@pytest.mark.parametrize("call, expected", LAMBERT)
def test_lambert_solutions(call, expected):
    assert sorted(str(s) for s in call()) == sorted(expected)


def test_lambert_solutions_check_numerically():
    for eq, sols in [(x*exp(x**2) - 1, solve(x*exp(x**2) - 1, x)),
                     (exp(x) + x**3, solve(exp(x) + x**3, x)),
                     (x**2 - 2**x, solve(x**2 - 2**x, x))]:
        for s in sols:
            assert abs(complex(N(eq.subs(x, s), 30))) < 1e-20


SOLVESET = [
    (lambda: solveset(sqrt(x)*sqrt(x - 7) - 12, x), '{16}'),
    (lambda: solveset(sqrt(x)*sqrt(x - 7) - 12, x, S.Reals), '{16}'),
    (lambda: solveset(2*sqrt(x)*sqrt(x - 7) - 24, x), '{16}'),
    (lambda: solveset(sqrt(x)*sqrt(x + 1) + 2, x), '{-sqrt(17)/2 - 1/2}'),
    (lambda: solveset(1/(1 + a*x) - 1, x), 'Complement({0}, {-1/a})'),
    (lambda: solveset(a/(x - 1) - 2, x), 'Complement({a/2 + 1}, {1})'),
    (lambda: solveset(x/(x - a), x), 'Complement({0}, {a})'),
    (lambda: solveset((x**2 - a)/(x - b), x),
     'Complement({-sqrt(a), sqrt(a)}, {b})'),
    (lambda: solveset((x - a)/(x - b), x), 'Complement({a}, {b})'),
    (lambda: solveset(1/(x - 1) - 2, x), '{3/2}'),
    (lambda: solveset(1/(x - a) - 1/(x - b), x), 'EmptySet'),
    (lambda: solveset(x**3 + x + 1, x, S.Reals),
     '{-(27/2 + 3*sqrt(93)/2)**(1/3)/3 + (27/2 + 3*sqrt(93)/2)**(-1/3)}'),
]


@pytest.mark.parametrize("call, expected", SOLVESET)
def test_solveset_forms(call, expected):
    assert str(call()) == expected


def test_solveset_symbolic_pole_complement():
    # SymPy: Complement({(b/y - 1)/a}, {-1/a}); the solution is written
    # b/(a*y) - 1/a by symox (same value)
    r = solveset(y - b/(1 + a*x), x)
    assert type(r).__name__ == 'Complement'
    sol, = list(r.args[0])
    assert (sol - (b/y - 1)/a).expand() == 0
    assert str(r.args[1]) == '{-1/a}'


def test_solveset_symbolic_pole_reals():
    # SymPy 1.14 gives the same Complement forms here
    assert str(solveset((x - a)/(x - b), x, S.Reals)) == \
        'Complement(Intersection({a}, Reals), Intersection({b}, Reals))'
    assert str(solveset(x/(x - a), x, S.Reals)) == \
        'Complement({0}, Intersection({a}, Reals))'
    assert str(solveset((x**2 - 1)/(x - a), x, S.Reals)) == \
        'Complement({-1, 1}, Intersection({a}, Reals))'


def test_solveset_symbolic_pole_reals_intentional_complement():
    # SymPy 1.14 gives Intersection({(b/y - 1)/a}, Reals), which holds the
    # pole x = -1/a when b = 0; symox removes the poles as in Complexes
    r = solveset(y - b/(1 + a*x), x, S.Reals)
    assert type(r).__name__ == 'Complement'
    assert str(r.args[1]) == 'Intersection({-1/a}, Reals)'


def test_roots_integer_basis_form():
    c = '(114*sqrt(12657)/78125 + 12459439/52734375)**(1/3)'
    expected = {
        '-1484/375 + 172564/(140625*%s) + 4*%s' % (c, c),
        '-1484/375 + 172564/(140625*(-1/2 + sqrt(3)*I/2)*%s) + '
        '4*(-1/2 + sqrt(3)*I/2)*%s' % (c, c),
        '-1484/375 + 4*(-1/2 - sqrt(3)*I/2)*%s + '
        '172564/(140625*(-1/2 - sqrt(3)*I/2)*%s)' % (c, c),
    }
    p = 3125*x**3 + 37100*x**2 + 100800*x - 82944
    r = roots(p, x)
    assert set(r.values()) == {1}
    assert {str(k) for k in r} == expected
    for k in r:
        assert abs(complex(N(p.subs(x, k), 40))) < 1e-25


def test_limit_of_log_singular_difference_is_never_nan():
    # SymPy gives -log(2); symox may leave the limit unevaluated but must
    # not return nan
    r = limit(-Ei(-2*x) + Ei(-x), x, 0, '+')
    assert r != S.NaN
    if not r.has(Limit):
        assert abs(complex(N(r + log(2)))) < 1e-12


def test_reals_never_false_empty_for_transcendental_candidates():
    # SymPy: {3}; symox misses the real inverter of a negative base and
    # must not turn its complex candidate into EmptySet
    from symox import Eq
    r = solveset(Eq((-5)**(x + 1), 625), x, S.Reals)
    assert r != S.EmptySet
    assert str(r).startswith('Intersection(')


def test_bivariate_type_dependent_generators():
    # SymPy: []
    assert solve(x**R(1, 3) + sqrt(x) + 1, x) == []


def test_nonreal_inequality_still_raises():
    # SymPy: Interval(-oo, oo), but (1 + I)**3 = -2 + 2*I has no order
    import pytest as _pytest
    from symox import I, solve_univariate_inequality
    assert str(solveset((I*x + 1)**3 - 2, x, S.Reals)) == 'EmptySet'
    with _pytest.raises(NotImplementedError):
        solve_univariate_inequality((I*x + 1)**3 < 2, x, relational=False)
    assert str(solve_univariate_inequality((x + I)*(x - I) > 2, x,
                                           relational=False)) == \
        'Union(Interval.open(-oo, -1), Interval.open(1, oo))'
