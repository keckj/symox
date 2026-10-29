"""Native solveset result checks and the native modular step.

Expected values are SymPy 1.14 results (``str`` of the result, so the
bound ``_n`` of an ImageSet compares by name), except where noted as an
intentional correction of a wrong reference result.
"""
import pytest

from symox import (
    Eq, I, Interval, Mod, Range, Rational, S, Symbol, log, solveset, sqrt,
)
from symox.polys.polyerrors import NotInvertible

x = Symbol('x')
y = Symbol('y')
k = Symbol('k', integer=True)


@pytest.mark.parametrize('f, domain, expected', [
    # `_check`: a solution where a subexpression of the equation is infinite
    (Eq(1/(x**2 - 4) + x, 1/(x**2 - 4) + 2), S.Reals, 'EmptySet'),
    (Eq(1/x + x, 1/x), S.Reals, 'EmptySet'),
    (x*log(x), S.Reals, '{1}'),
    (x*log(x), S.Complexes, '{1}'),
    (x**2*log(x)/x, S.Reals, '{1}'),
    # `f.expand().is_zero`
    (-x**2 - 2*x + (x + 1)**2 - 1, S.Reals, 'Reals'),
    (-x**2 - 2*x + (x + 1)**2 - 1, S.Complexes, 'Complexes'),
    ((x + 1)**2 - x**2 - 2*x - 1, Interval(0, 1), 'Interval(0, 1)'),
    (x*(x + 1) - x**2 - x, S.Reals, 'Reals'),
    (Eq((x + 1)**2, x**2 + 2*x + 1), S.Complexes, 'Complexes'),
    # a numerator that vanishes: the domain minus the poles
    ((x**2 - 1)/(x - 1) - x - 1, S.Complexes, 'Complement(Complexes, {1})'),
    # the singularities of the relational are excluded
    (x + 1/x > -2 + 1/x, S.Reals, 'Union(Interval.open(-2, 0), Interval.open(0, oo))'),
    (x + 1/x > -2 + 1/x, S.Integers, 'Union(Range(-1, 0, 1), Range(1, oo, 1))'),
    (log(x) + 1/log(x) > -2 + 1/log(x), S.Integers, 'Range(2, oo, 1)'),
    (x + 1/x >= 1/x, S.Reals, 'Interval.open(0, oo)'),
    (log(x)/x <= 0, S.Integers, 'Range(1, 2, 1)'),
    (x/(x - 3) < 3, Range(0, 5), 'Range(0, 3, 1)'),
])
def test_result_checks(f, domain, expected):
    assert str(solveset(f, x, domain)) == expected


def test_logarithmic_intentional_correction():
    # SymPy 1.14 also returns -sqrt(10), where the sum of the principal
    # logarithms is 2*pi*I, not 0
    assert solveset(log(x - 3) + log(x + 3), x, S.Reals) == {sqrt(10)}


@pytest.mark.parametrize('f, domain, expected', [
    (Mod(x, 5) - 3, S.Integers, 'ImageSet(Lambda(_n, 5*_n + 3), Integers)'),
    (Eq(Mod(x, 5), 3), S.Integers, 'ImageSet(Lambda(_n, 5*_n + 3), Integers)'),
    (2*Mod(x, 3) - 2, S.Integers, 'ImageSet(Lambda(_n, 3*_n + 1), Integers)'),
    (2*Mod(x, 5) - 3, S.Integers, 'EmptySet'),
    (7 - Mod(x, 5), S.Integers, 'EmptySet'),
    (2 + Mod(x, -3), S.Integers, 'ImageSet(Lambda(_n, -3*_n - 2), Integers)'),
    (3 - Mod(5*x - 8, 7), S.Integers, 'ImageSet(Lambda(_n, 7*_n + 5), Integers)'),
    (Mod(2*x + 1, 9) - 4, S.Integers, 'ImageSet(Lambda(_n, 9*_n + 6), Integers)'),
    (Mod(x**2, 7) - 2, S.Integers,
     'Union(ImageSet(Lambda(_n, 7*_n + 3), Integers), ImageSet(Lambda(_n, 7*_n + 4), Integers))'),
    (Mod(x**2, 9), S.Integers, 'Union(ImageSet(Lambda(_n, 9*_n), Integers), '
     'ImageSet(Lambda(_n, 9*_n + 3), Integers), ImageSet(Lambda(_n, 9*_n + 6), Integers))'),
    (Mod(x**3, 7) - 2, S.Integers, 'EmptySet'),
    (Mod(x**4, 14) - 11, S.Integers,
     'Union(ImageSet(Lambda(_n, 14*_n + 3), Integers), ImageSet(Lambda(_n, 14*_n + 11), Integers))'),
    (Mod(7**x, 41) - 15, S.Integers, 'ImageSet(Lambda(_n, 40*_n + 3), Naturals0)'),
    (Mod(12**x, 21) - 18, S.Integers, 'ImageSet(Lambda(_n, 6*_n + 2), Naturals0)'),
    (Mod((-2)**x, 7) - 1, S.Integers, 'ImageSet(Lambda(_n, 6*_n), Naturals0)'),
    (Mod(5**x, 25), S.Integers, 'ImageSet(Lambda(_n, 20*_n), Naturals0)'),
    (Mod(x, 3) - I, S.Integers, 'EmptySet'),
    (Mod(I*x, 3) - 2, S.Integers, 'ConditionSet(x, Eq(Mod(I*x, 3) - 2, 0), Integers)'),
    (Mod(sqrt(x), 5) - 1, S.Integers, 'ConditionSet(x, Eq(Mod(sqrt(x), 5) - 1, 0), Integers)'),
    (x - Mod(x, 5), S.Integers, 'ConditionSet(x, Eq(x - Mod(x, 5), 0), Integers)'),
    (Mod(x, 5) - 3, S.Reals, 'ConditionSet(x, Eq(Mod(x, 5) - 3, 0), Reals)'),
    (Mod(x, 5) - 3, S.Naturals, 'Intersection(ImageSet(Lambda(_n, 5*_n + 3), Integers), Naturals)'),
    (3 - Mod(5*x - 8, 7), S.Naturals0, 'Intersection(ImageSet(Lambda(_n, 7*_n + 5), Integers), Naturals0)'),
])
def test_modular(f, domain, expected):
    assert str(solveset(f, x, domain)) == expected


@pytest.mark.parametrize('f, error', [
    (Mod(x, k) - 1, TypeError),
    (Mod(x, 5) - y, TypeError),
    (Mod(x/2, 5) - 1, NotInvertible),
    (x*Mod(x, 5) - 1, AttributeError),
    (Mod(y**x, 5) - 1, AttributeError),
])
def test_modular_errors(f, error):
    with pytest.raises(error):
        solveset(f, x, S.Integers)


def test_modular_intentional_corrections():
    # SymPy 1.14 gives ImageSet(Lambda(_n, 5*_n), Integers), but k = x = 1
    # solves Mod(k*x, 5) = 1
    assert str(solveset(Mod(k*x, 5) - 1, x, S.Integers)) == \
        'ConditionSet(x, Eq(Mod(k*x, 5) - 1, 0), Integers)'
    # SymPy 1.14 uses the numerator of the base: 4*Naturals0, but
    # Mod((1/2)**4, 5) = 1/16
    assert str(solveset(Mod(Rational(1, 2)**x, 5) - 1, x, S.Integers)).startswith('ConditionSet(')
