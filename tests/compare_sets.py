"""Differential tests for sets, ``solveset`` and inequality solving against
SymPy.

Each case is evaluated with SymPy and symox; printed forms are compared,
and when they differ results are compared after parsing the symox output
with SymPy (set equality / logical equivalence).
"""

from __future__ import annotations

import sys

import sympy
import symox

CASES = r"""
Interval(0, 1)
Interval(0, 1, True)
Interval(0, 1, False, True)
Interval.open(0, 1)
Interval.Lopen(0, 1)
Interval.Ropen(0, 1)
Interval(1, 0)
Interval(1, 1)
Interval(1, 1, True)
Interval(-oo, oo)
Interval(-oo, 0)
Interval(0, oo)
Interval(-oo, 0, True)
Interval(x, y)
Interval(0, x)
FiniteSet(3, 1, 2)
FiniteSet(x, 1, y)
FiniteSet(y, x)
FiniteSet(1, 1, 2)
FiniteSet()
FiniteSet(pi, 1, Rational(1,2))
FiniteSet(Interval(0,1))
FiniteSet(1, x, Rational(1,2), -1, pi, sqrt(2), I)
Union(Interval(0, 1), Interval(2, 3))
Union(Interval(2, 3), Interval(0, 1))
Union(Interval(0, 2), Interval(1, 3))
Union(Interval(0, 1), Interval(1, 2))
Union(Interval(0, 1, False, True), Interval(1, 2, True))
Union(Interval(0, 1, False, True), Interval(1, 2))
Union(Interval(0, 1), FiniteSet(2))
Union(Interval(0, 1), FiniteSet(1))
Union(Interval(0, 1, True, True), FiniteSet(0, 1))
Union(Interval.open(0, 1), FiniteSet(0))
Union(FiniteSet(1,2), FiniteSet(3))
Union(FiniteSet(x), FiniteSet(1))
Union(Interval(0,1), FiniteSet(x))
Union(Interval(0,1), FiniteSet(5, x))
Intersection(Interval(0, 2), Interval(1, 3))
Intersection(Interval(0, 1), Interval(2, 3))
Intersection(Interval(0, 1), Interval(1, 2))
Intersection(Interval(0, 1), FiniteSet(0, 2, Rational(1,2)))
Intersection(Interval(0, 1), FiniteSet(x))
Intersection(FiniteSet(1,2,3), FiniteSet(2,3,4))
Intersection(Union(Interval(0,1), Interval(2,3)), Interval(Rational(1,2), Rational(5,2)))
Complement(Interval(0, 2), Interval(1, 3))
Complement(Interval(0, 3), Interval(1, 2))
Complement(Interval(0, 3), FiniteSet(1))
Complement(FiniteSet(1, 2, 3), FiniteSet(2))
Complement(S.Reals, FiniteSet(0))
Complement(S.Reals, Interval(0, 1))
Complement(FiniteSet(x, 1), FiniteSet(1))
Complement(FiniteSet(x, 1), Interval(0, 2))
Complement(Interval(0,1), FiniteSet(x))
S.EmptySet
S.UniversalSet
S.Reals
S.Integers
S.Naturals
S.Naturals0
S.Complexes
S.Rationals
Union(S.EmptySet, Interval(0,1))
Intersection(S.Reals, Interval(0,1))
Intersection(S.Integers, Interval(0, 5))
Intersection(S.Integers, Interval(Rational(1,2), Rational(7,2)))
Intersection(S.Naturals, Interval(0, 3))
Intersection(S.Integers, Interval(0, oo))
Union(S.Reals, Interval(0,1))
Interval(0, 1).inf
Interval(0, 1).sup
Interval(0, 1).measure
Union(Interval(0, 1), Interval(2, 4)).measure
FiniteSet(1,2).measure
Interval(0, 1).boundary
Interval(0, 1, True).boundary
FiniteSet(1, 2).boundary
Union(Interval(0, 1), Interval(2, 3)).boundary
Interval(0, oo).boundary
Interval(0, 1).contains(x)
Interval(0, 1, True, True).contains(x)
Interval(0, 1).contains(Rational(1,2))
Interval(0, 1).contains(2)
FiniteSet(1, 2).contains(x)
FiniteSet(1, 2).contains(1)
Union(Interval(0, 1), Interval(2, 3)).contains(x)
Interval(0, oo).contains(x)
S.Reals.contains(x)
S.Integers.contains(x)
S.Integers.contains(2)
S.Integers.contains(Rational(1,2))
Interval(0, 1).as_relational(x)
Interval(0, 1, True, True).as_relational(x)
Interval(0, oo).as_relational(x)
Interval(-oo, 1).as_relational(x)
Interval(-oo, oo).as_relational(x)
FiniteSet(1, 2).as_relational(x)
Union(Interval(0, 1), Interval(2, 3)).as_relational(x)
Union(Interval(0, 1), FiniteSet(3)).as_relational(x)
Complement(S.Reals, FiniteSet(0)).as_relational(x)
Interval(0, 1) + Interval(2, 3)
Interval(0, 1) | Interval(2, 3)
Interval(0, 2) & Interval(1, 3)
Interval(0, 2) - Interval(1, 3)
FiniteSet(1,2) - FiniteSet(1)
Interval(0,1).is_subset(Interval(0,2))
Interval(0,1).union(Interval(0.5,2))
Interval(0,1).intersect(Interval(0.5,2))
Interval(0,1).complement(S.Reals)
Interval(0, 1).is_open
Interval(0,1).left_open
Interval(0, 1).left
Interval(0, 1).end
ProductSet(Interval(0,1), FiniteSet(1,2))
Interval(0,1)**2
ImageSet(Lambda(n, 2*n*pi), S.Integers)
ImageSet(Lambda(n, n**2), S.Naturals)
ConditionSet(x, Eq(sin(x), x), S.Reals)
Contains(x, Interval(0,1))
Interval(0,1).contains(x).func
FiniteSet(1, 2, 3).inf
FiniteSet(x, 2).inf
Union(Interval(0,1), Interval(3,4)).sup
S.Reals.inf
S.EmptySet.measure
Interval(x, x+1).measure
Interval(0, 1) == Interval(0, 1)
2 in Interval(0,3)
Interval(0, 1).is_empty
FiniteSet(1, 2, 3).is_empty
Union(FiniteSet(1, 2), Interval(3, 4))
Union(FiniteSet(1, 2), Interval(1, 4))
Intersection(Interval(0, 2), Interval(1, 3), Interval(Rational(3,2), 4))
Union(Interval(0, 1), Interval(2, 3), Interval(Rational(1,2), Rational(5,2)))
Complement(Interval(0, 10), Union(Interval(1, 2), Interval(3, 4)))
Interval(0, 1, True, True).closure
Interval(0, 1).interior
Intersection(FiniteSet(x, y), FiniteSet(x))
Intersection(Interval(0, 1), Interval(x, 2))
Union(Interval(x, 1), Interval(2, 3))
Range(0, 0)
Range(1, 2)
Range(0, 6)
Range(5)
Range(0, 10, 2)
Range(0, 6).inf
Range(0, 6).sup
Range(0, 6).contains(x)
Range(0, 6).contains(3)
Range(0, 6).measure
Intersection(S.Integers, Interval(-oo, 2))
Intersection(S.Integers, Interval(Rational(1,2), 1))
Intersection(S.Integers, Interval(Rational(1,2), Rational(3,4)))
Intersection(S.Integers, Interval.open(0, 3))
Intersection(S.Naturals0, Interval(-2, 3))
Intersection(S.Integers, Interval(x, 3))
S.Complexes.contains(x)
S.Complexes.contains(1)
S.Reals.contains(xr)
S.Reals.contains(I)
S.Integers.contains(Symbol('n', integer=True))
S.Naturals.contains(0)
S.Naturals0.contains(0)
S.Rationals.contains(Rational(1,2))
S.Rationals.contains(pi)
S.EmptySet.inf
S.Integers.inf
S.Naturals.inf
S.Naturals.sup
S.Reals.measure
S.Integers.measure
S.Reals.boundary
S.Integers.boundary
Intersection(FiniteSet(x, y), FiniteSet(x, z))
Intersection(FiniteSet(1, x), FiniteSet(1, 2))
Intersection(FiniteSet(1, x), FiniteSet(2, 3))
Intersection(FiniteSet(1, 2, x), Interval(0, Rational(3,2)))
Intersection(FiniteSet(x, y), Interval(0, 1))
Union(FiniteSet(1), Interval.open(0, 1))
Union(FiniteSet(0, 1, 5), Interval.open(0, 1))
Union(Interval(0, 1), Interval(x, 2))
Complement(FiniteSet(x, 1, 2), FiniteSet(1, y))
Complement(FiniteSet(1, 2, 3), Interval(0, Rational(5,2)))
Complement(Interval(0, 1), Interval(0, 1))
Complement(Interval(0, 1), S.Reals)
Complement(S.Reals, S.Integers)
Complement(S.Integers, FiniteSet(1))
Complement(Interval(0, 5), FiniteSet(1, x))
Complement(S.Reals, FiniteSet(1, 2, x))
Complement(Union(Interval(0,1), Interval(2,3)), FiniteSet(Rational(1,2)))
Intersection(S.Reals, S.Integers)
Intersection(S.Complexes, Interval(0,1))
Union(S.Integers, S.Reals)
Union(S.Integers, FiniteSet(Rational(1,2)))
Union(S.Integers, FiniteSet(2))
Intersection(Interval(0,1), Interval(0,1,True,True))
Intersection(Interval(0,2,True), Interval(0,2,False,True))
Union(Interval(0,1,True,True), Interval(1,2,True,True))
Union(Interval(0,1,True,True), Interval(1,2,True,True), FiniteSet(1))
Interval(0, 1).is_subset(Interval(0, 1, True))
Interval(0, 1).is_subset(S.Reals)
FiniteSet(1, 2).is_subset(Interval(0, 3))
FiniteSet(x).is_subset(Interval(0, 3))
Interval(0, 2).is_superset(FiniteSet(1))
Interval(0, 2).is_proper_subset(Interval(0, 3))
Interval(0, oo).measure
Interval(0, 1).is_left_unbounded
Interval(0, 1).is_finite_set
FiniteSet(1,2,3).is_finite_set
len(FiniteSet(1,2,3))
list(FiniteSet(3,2,1))
Union(Interval(0,1), Interval(2,3)).is_iterable
ImageSet(Lambda(x, x**2), FiniteSet(1,2,3))
ImageSet(Lambda(x, x), S.Integers)
ConditionSet(x, x > 0, FiniteSet(-1, 1, y))
ConditionSet(x, x > 0, S.Reals)
ConditionSet(x, Eq(x**2, 1), S.UniversalSet)
ProductSet(Interval(0,1), S.EmptySet)
ProductSet(Interval(0,1))
Interval(0, 1).complement(Interval(0, 2))
FiniteSet(1, 2).complement(S.Reals)
FiniteSet(1,2) | Interval(0, 1)
Interval(0, 1) ^ Interval(0, 2)
SymmetricDifference(Interval(0, 2), Interval(1, 3))
Interval(0, 1).union(FiniteSet(2)).inf
Interval(1, 0, True)
Interval(oo, oo)
Interval(-oo, -oo)
Interval(0, 1).closure
Union(Interval(0, 1), Interval(2, 3)).closure
FiniteSet(1, 2).closure
Union(Interval.Ropen(0, 1), Interval.Lopen(1, 2)).boundary
Union(Interval.open(0, 1), FiniteSet(3)).boundary
Union(Interval(0,1), FiniteSet(3)).measure
Interval(0, 1).contains(I)
Interval(0, 1).contains(oo)
Interval(0, oo).contains(oo)
Contains(2, Interval(0, 3))
Contains(x, S.Reals)
Contains(I, S.Reals)
solveset(x**2 - 1, x)
solveset(x**2 - 1, x, S.Reals)
solveset(x**2 + 1, x)
solveset(x**2 + 1, x, S.Reals)
solveset(x - 1, x)
solveset(x**3 - x, x)
solveset(x**3 - 1, x)
solveset(x**3 - 1, x, S.Reals)
solveset(1/x, x)
solveset(x/(x-1), x)
solveset((x-1)/(x-1), x)
solveset(exp(x) - 1, x)
solveset(exp(x) - 1, x, S.Reals)
solveset(exp(x) - 2, x, S.Reals)
solveset(exp(x) - 2, x)
solveset(sin(x), x)
solveset(sin(x), x, S.Reals)
solveset(sin(x) - 1, x, S.Reals)
solveset(cos(x), x, S.Reals)
solveset(cos(x) - 1, x, S.Reals)
solveset(tan(x), x, S.Reals)
solveset(sin(x) - Rational(1,2), x, S.Reals)
solveset(sin(x) - x, x, S.Reals)
solveset(x - x, x)
solveset(x - x + 1, x)
solveset(Eq(x, 2), x)
solveset(x > 1, x, S.Reals)
solveset(x**2 > 1, x, S.Reals)
solveset(x**2 < 4, x, S.Reals)
solveset(x**2 <= 0, x, S.Reals)
solveset(log(x) - 1, x, S.Reals)
solveset(x**2 - 2*x + 1, x)
solveset(x**4 - 1, x)
solveset(x**4 - 1, x, S.Reals)
solveset(x**2 - y, x)
solveset(Abs(x) - 1, x, S.Reals)
solveset(sqrt(x) - 2, x)
solveset(x**2 - 1, x, Interval(0, oo))
solve_univariate_inequality(x**2 > 1, x)
solve_univariate_inequality(x**2 > 1, x, relational=False)
solve_univariate_inequality(x**2 - 3*x + 2 < 0, x)
solve_univariate_inequality(x**2 - 3*x + 2 <= 0, x, relational=False)
solve_univariate_inequality(x**2 + 1 > 0, x)
solve_univariate_inequality(x**2 + 1 < 0, x)
solve_univariate_inequality(x**3 - x >= 0, x, relational=False)
solve_univariate_inequality((x-1)/(x+1) > 0, x, relational=False)
solve_univariate_inequality((x-1)/(x+1) >= 0, x)
solve_univariate_inequality(1/x > 0, x)
solve_univariate_inequality(x > 0, x)
solve_univariate_inequality(x >= 0, x)
solve_univariate_inequality(x < 0, x)
solve_univariate_inequality(2*x - 3 > 0, x)
solve_univariate_inequality(x**2 >= 0, x)
solve_univariate_inequality(x**2 > 0, x)
solve_univariate_inequality(x**2 < 0, x)
solve_univariate_inequality(Ne(x**2, 1), x)
solve_univariate_inequality(x**2 - 2 > 0, x)
solve_univariate_inequality((x-2)**2*(x+1) > 0, x)
reduce_inequalities(x**2 > 1, x)
reduce_inequalities([x > 1, x < 3], x)
reduce_inequalities(x**2 - 3*x + 2 < 0, x)
reduce_inequalities([x**2 > 1, x < 3], x)
reduce_inequalities(x > 1)
reduce_inequalities(2*x + 1 >= 0, x)
reduce_inequalities((x-1)/(x-2) <= 0, x)
reduce_inequalities(x**2 < 0, x)
reduce_inequalities(x**2 >= 0, x)
reduce_inequalities(Eq(x**2, 1), x)
reduce_inequalities(r**2 > 1)
solveset(sin(2*x), x, S.Reals)
solveset(sin(x + 1), x, S.Reals)
solveset(cos(x/2), x, S.Reals)
solveset(sin(x) + 1, x, S.Reals)
solveset(tan(2*x) - 1, x, S.Reals)
solveset(sin(x) - sqrt(2)/2, x, S.Reals)
solveset(2*sin(x) - 1, x, S.Reals)
solveset(exp(2*x) - 1, x)
solveset(exp(x) + 1, x)
solveset(exp(x) + 1, x, S.Reals)
solveset(cos(x) + 1, x, S.Reals)
solveset(cos(x) - Rational(1,2), x, S.Reals)
solveset(sin(x) - 2, x, S.Reals)
solveset(sin(x)*cos(x), x, S.Reals)
solveset(sin(x) - y, x, S.Reals)
solveset(tan(x) - 1, x, S.Reals)
solveset(cos(2*x + 1), x, S.Reals)
solveset(sin(x), x, Interval(0, 10))
solveset(y, x)
solveset(x**2 - 4, x, Interval(0, 5))
solveset(x**2 + y, x, S.Reals)
solveset(x*(x-1), x, S.Reals)
solveset(x**5 - x - 1, x, S.Reals)
solveset(log(x), x)
solveset(exp(x), x)
solveset(Ne(x, 1), x, S.Reals)
solveset(x**2 > 4, x)
solveset(Abs(x) - 2, x)
reduce_inequalities(x**2 + 1 > 0, x)
reduce_inequalities(x**2 > 0, x)
reduce_inequalities(Ne(x, 1), x)
reduce_inequalities([x > 1, x > 2], x)
reduce_inequalities([x > 1, x < 0], x)
reduce_inequalities(x < 2, x)
reduce_inequalities(x**3 - x > 0, x)
reduce_inequalities(1/x > 0, x)
reduce_inequalities((x+1)/(x-1) >= 0, x)
reduce_inequalities([x >= 0, x <= 5, Ne(x, 3)], x)
reduce_inequalities(x**2 <= 0, x)
solve_univariate_inequality(x**2 <= 0, x)
solve_univariate_inequality(Eq(x, 1), x)
solve_univariate_inequality(x**2 - 4 <= 0, x, relational=False, domain=Interval(0, 10))
solve_univariate_inequality(exp(x) > 1, x)
solve_univariate_inequality(sqrt(x) < 1, x)
reduce_inequalities(Abs(x) < 1, x)
reduce_inequalities(x > y, x)
Interval(0, 1).contains(Rational(1, 3))
Interval(0, pi).contains(3)
Interval(0, pi).contains(4)
Interval(sqrt(2), 2).contains(Rational(3, 2))
Interval(0, 1) - FiniteSet(Rational(1, 2))
Interval(-oo, oo) - Interval(0, oo)
Union(Interval(0, 1), Interval(Rational(1, 2), 2), FiniteSet(5))
Intersection(Interval(0, 5), Union(FiniteSet(1, 7), Interval(3, 9)))
Complement(Interval(0, 10), Interval(2, 3), evaluate=True) if False else Complement(Interval(0, 10), Interval(2, 3))
FiniteSet(1, 2, 3) - Interval(2, 3)
FiniteSet(1, 2, 3) & Interval(2, 3)
FiniteSet(1, 2, 3) | FiniteSet(3, 4)
FiniteSet(x, y) - FiniteSet(x)
FiniteSet(Rational(1, 3), Rational(1, 2), 0.25)
FiniteSet(-oo, oo, 0)
FiniteSet(FiniteSet(1, 2), 3)
Interval(0, 1) * FiniteSet(1, 2)
Interval(0, 2).measure
Interval(0, 2, True, True).measure
Union(Interval(0, 1), Interval(3, 5)).measure
Union(Interval(0, 1), Interval(3, 5)).inf
Union(Interval(0, 1), Interval(3, 5)).sup
FiniteSet(3, 1, 2).sup
Interval(0, 1).is_subset(Interval(-1, 2))
Interval(0, 2).is_subset(Interval(0, 1))
Interval(0, 1).is_subset(FiniteSet(0, 1))
FiniteSet(0, 1).is_subset(Interval(0, 1))
S.Integers.is_subset(S.Reals)
S.Reals.is_subset(S.Integers)
S.Naturals.is_subset(S.Naturals0)
Union(Interval(0, 1), Interval(2, 3)).is_subset(Interval(0, 3))
Interval(0, 1).union(Interval(2, 3)).as_relational(x)
Interval.Ropen(0, 1).as_relational(x)
Interval.Lopen(-oo, 1).as_relational(x)
Union(Interval.Ropen(0, 1), Interval.Lopen(1, 2)).as_relational(x)
Union(Interval.open(0, 1), Interval.open(1, 2), Interval.open(2, 3)).as_relational(x)
Intersection(Interval(0, 1), FiniteSet(x)).as_relational(x)
Interval(0, 1).boundary
Interval(0, oo).closure
S.Reals.contains(1)
S.Reals.contains(oo)
S.Naturals0.contains(-1)
S.Integers.contains(pi)
S.Rationals.contains(sqrt(2))
solveset(x**2 - 5*x + 6, x)
solveset(x**2 + 2*x + 5, x)
solveset(x**2 + 2*x + 5, x, S.Reals)
solveset(x**4 - 5*x**2 + 4, x, S.Reals)
solveset(x**3 + x, x)
solveset(x**3 + x, x, S.Reals)
solveset((x - 1)*(x - 2)/(x - 1), x)
solveset(x**2/x, x)
solveset(exp(x) - 3, x, S.Reals)
solveset(log(x) - 2, x)
solveset(2**x - 8, x, S.Reals)
solveset(sqrt(x) - 3, x, S.Reals)
solveset(cos(x) - 2, x, S.Reals)
solveset(tan(x) - sqrt(3), x, S.Reals)
solveset(sin(3*x) - 1, x, S.Reals)
solveset(cos(x) + Rational(1, 2), x, S.Reals)
solveset(2*cos(x)**2 - 1, x, S.Reals)
solveset(sin(x)**2 - 1, x, S.Reals)
solveset(x - 3 > 0, x, S.Reals)
solveset(x**2 - 1 <= 0, x, S.Reals)
solveset((x - 1)/(x + 2) < 0, x, S.Reals)
solveset(Eq(x**2, 4), x, S.Reals)
solveset(Ne(x**2, 4), x, S.Reals)
solveset(x**2 - 1, x, Interval(-5, 0))
solveset(x**2 - 1, x, Interval.open(-1, 1))
solveset(cos(x), x, Interval(0, 2*pi))
solveset(sin(2*x), x, Interval(0, pi))
solve_univariate_inequality(x**2 - 1 < 0, x, relational=False)
solve_univariate_inequality(x**3 - 3*x + 1 > 0, x, relational=False)
solve_univariate_inequality(x**4 - 1 >= 0, x)
solve_univariate_inequality((x - 1)**2 >= 0, x)
solve_univariate_inequality((x - 1)**2 > 0, x)
solve_univariate_inequality((x - 1)**2 <= 0, x)
solve_univariate_inequality(x/(x**2 + 1) > 0, x)
solve_univariate_inequality((x**2 - 4)/(x - 1) <= 0, x, relational=False)
solve_univariate_inequality(1/(x - 1) < 1, x, relational=False)
solve_univariate_inequality(x**2 - 2*x - 1 > 0, x, relational=False)
solve_univariate_inequality(log(x) > 0, x)
solve_univariate_inequality(exp(x) < 2, x, relational=False)
solve_univariate_inequality(sqrt(x) > 2, x, relational=False)
solve_univariate_inequality(Abs(x) < 2, x, relational=False)
solve_univariate_inequality(Abs(x - 1) >= 3, x, relational=False)
solve_univariate_inequality(x**2 > 1, x, relational=False, domain=Interval(0, 10))
reduce_inequalities(x**2 - 4 > 0, x)
reduce_inequalities([x**2 - 4 > 0, x > 0], x)
reduce_inequalities([x**2 < 9, x > -1], x)
reduce_inequalities([x >= 1, x <= 1], x)
reduce_inequalities(Abs(x - 2) <= 1, x)
reduce_inequalities(x**3 > 8, x)
reduce_inequalities(2*x - 4 < x, x)
reduce_inequalities(-x > 2, x)
reduce_inequalities(x**2 + x + 1 > 0, x)
reduce_inequalities(x**2 + x + 1 < 0, x)
""".strip().splitlines()

# Deliberate deviations (documented in the report).
KNOWN_DEVIATIONS = {
    # SymPy's CRootOf is not implemented in symox
    "solveset(x**5 - x - 1, x, S.Reals)",
    # SymPy raises; symox returns a ConditionSet
    "solveset(Abs(x) - 2, x)",
    # SymPy factors I out of the ImageSet body
    "solveset(exp(x) + 1, x)",
    # SymPy's printing order depends on Dummy creation order
    "solveset(sin(x)*cos(x), x, S.Reals)",
    # SymPy raises TypeError; symox keeps the intersection unevaluated
    "Intersection(S.Integers, Interval(x, 3))",
}


def ns(mod):
    d = {n: getattr(mod, n) for n in dir(mod) if not n.startswith("_")}
    if mod is sympy:
        from sympy.solvers.inequalities import reduce_inequalities, solve_univariate_inequality
        d["reduce_inequalities"] = reduce_inequalities
        d["solve_univariate_inequality"] = solve_univariate_inequality
    d.update(zip("xyzn", mod.symbols("x y z n")))
    d["r"] = mod.Symbol("x", real=True)
    return d


def semantic_equal(got, want):
    try:
        g = sympy.sympify(str(got))
    except Exception:  # noqa: BLE001
        return False
    if g == want:
        return True
    try:
        if isinstance(want, sympy.Set) and isinstance(g, sympy.Set):
            return bool(sympy.Complement(g, want) == sympy.S.EmptySet and sympy.Complement(want, g) == sympy.S.EmptySet)
        if isinstance(want, sympy.logic.boolalg.Boolean):
            return sympy.simplify(sympy.Equivalent(g, want)) == sympy.true
    except Exception:  # noqa: BLE001
        return False
    return False


def compare(case, sx, mx):
    try:
        want = eval(case, sx)
        want_err = None
    except Exception as e:  # noqa: BLE001
        want, want_err = None, type(e).__name__
    try:
        got = eval(case, mx)
        got_err = None
    except Exception as e:  # noqa: BLE001
        got, got_err = None, f"{type(e).__name__}: {e}"
    if want_err is not None:
        if got_err is not None:
            return "both-error", got_err, want_err
        return "DIFFERENT", str(got), f"error {want_err}"
    if got_err is not None:
        return "error", got_err, str(want)
    if str(got) == str(want):
        return "same", str(got), str(want)
    if semantic_equal(got, want):
        return "equivalent", str(got), str(want)
    return "DIFFERENT", str(got), str(want)


def run(verbose=True):
    sx, mx = ns(sympy), ns(symox)
    bad, stats = [], {}
    for c in CASES:
        st, got, want = compare(c, sx, mx)
        stats[st] = stats.get(st, 0) + 1
        if st not in ("same", "both-error"):
            if verbose:
                print(f"[{st:10s}] {c:55s} symox: {got:40s} sympy: {want}")
            if st in ("DIFFERENT", "error") and c not in KNOWN_DEVIATIONS:
                bad.append((c, got, want))
    if verbose:
        print("\nsummary:", stats, f"(total {len(CASES)})")
    return bad


def test_sets_against_sympy():
    bad = run(verbose=False)
    assert not bad, "\n".join(f"{c}: symox={g} sympy={w}" for c, g, w in bad)


if __name__ == "__main__":
    sys.exit(1 if run() else 0)
