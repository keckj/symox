"""Union emptiness and the subset decision path match the reference.

Expected values were produced by the reference SymPy 1.14 on the same
expressions (Union.is_empty is the fuzzy and of its components,
is_subset follows trivial cases, the class hooks, the pairwise rules and
the intersection fallback).
"""
import pytest

import symox as sx
from symox import (S, Interval, Union, Intersection, Complement, FiniteSet,
                   ProductSet, Range, pi, oo, sin, is_increasing, is_decreasing)

x, y, a, b = sx.symbols('x y a b', real=True)

U = Union(Interval(0, 1), Interval(2, 3))
W = Union(Interval(0, pi/2), Interval.Ropen(3*pi/2, 2*pi))

CASES = [
    ('U.is_empty', False),
    ('U.is_finite_set', False),
    ('Union(FiniteSet(1), Interval(2,3)).is_empty', False),
    ('Union(FiniteSet(x), Interval(a,b)).is_empty', False),
    ('Interval(a,b).is_empty', None),
    ('Interval(a,a+1).is_empty', False),
    ('Interval(a,b).is_finite_set', None),
    ('Intersection(Interval(0,1), Interval(a,b)).is_empty', None),
    ('Complement(Interval(0,1), Interval(a,b)).is_empty', None),
    ('ProductSet(Interval(0,1), Interval(2,3)).is_empty', False),
    ('Intersection(Interval(0,1), Interval(a,b)).is_finite_set', None),
    ('Complement(Interval(0,1), FiniteSet(a)).is_finite_set', False),
    ('S.Reals.is_subset(W)', False),
    ('Interval(0,10).is_subset(W)', False),
    ('Interval(0,1).is_subset(W)', True),
    ('Interval(1,2).is_subset(W)', None),
    ('U.is_subset(Interval(0,3))', True),
    ('U.is_subset(Interval(0,2))', False),
    ('U.is_subset(S.Reals)', True),
    ('Interval(0,3).is_superset(U)', True),
    ('Interval(0,2).is_superset(U)', False),
    ('U.is_superset(Interval(0,1))', True),
    ('U.is_superset(Interval(0,3))', None),
    ('U.is_proper_subset(U)', False),
    ('U.is_proper_subset(Interval(0,3))', True),
    ('Interval(0,1).is_proper_subset(U)', True),
    ('U.is_proper_superset(FiniteSet(0,3))', True),
    ('U.is_subset(U)', True),
    ('Interval(0,1).is_subset(Union(Interval(0,S(1)/2), Interval(S(1)/2,1)))', True),
    ('Interval(a,b).is_subset(FiniteSet(1))', None),
    ('Interval(0,1).is_subset(FiniteSet(1,2))', False),
    ('Interval(0,1).is_subset(Interval(a,b))', None),
    ('FiniteSet(1,2).is_subset(Union(Interval(0,1),FiniteSet(2)))', True),
    ('FiniteSet(x).is_subset(Interval(0,1))', None),
    ('S.Naturals.is_subset(Interval(0,oo))', True),
    ('S.Integers.is_subset(S.Reals)', True),
    ('S.Rationals.is_subset(Interval(0,1))', None),
    ('Range(1,5).is_subset(Interval(0,4))', True),
    ('Range(1,5).is_subset(Interval.open(1,5))', False),
    ('Range(0,3).is_subset(FiniteSet(0,1,2))', True),
    ('Interval(0,1).is_subset(S.Rationals)', False),
    ('S.Complexes.is_subset(S.Reals)', False),
    ('Union(Interval(0,1),FiniteSet(5)).is_subset(Interval(0,4))', False),
    ('Intersection(Interval(0,1),Interval(a,b)).is_subset(Interval(0,1))', True),
    ('is_increasing(sin(x), symbol=x)', False),
    ('is_decreasing(sin(x), symbol=x)', False),
    ('is_increasing(x, Union(Interval(0,1),Interval(2,3)), x)', True),
    ('is_increasing(x, FiniteSet(0,1), x)', True),
    ('is_increasing(x, S.Complexes, x)', False),
]


@pytest.mark.parametrize("expression,expected", CASES)
def test_reference_values(expression, expected):
    assert eval(expression) is expected


def test_superset_requires_a_set():
    with pytest.raises(ValueError):
        U.is_superset(1)
    with pytest.raises(ValueError):
        Interval(0, 1).is_subset(1)
    with pytest.raises(ValueError):
        U.is_proper_subset(1)


def test_whole_line_intersection_keeps_incomparable_interval():
    z = sx.Symbol('z')
    assert Interval(0, z).intersect(S.Reals) == Interval(0, z)
    assert Interval(0, z).is_subset(S.Reals) is True
