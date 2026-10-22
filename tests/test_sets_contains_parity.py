"""Reference membership, inversion and diophantine rules of SymPy 1.14:
`Range.contains` with any step, `ImageSet.contains` and `ProductSet.contains`
of symbolic elements, `solveset`/`invert_complex`/`invert_real` ignoring the
assumptions of the symbol, the parabolic binary quadratic and the reference
linear parameterization of `diophantine`, and the image-set intersections
that rely on them.

Expected values are the srepr of the reference SymPy 1.14 on the same calls
(diophantine solutions as the sorted strings of their tuples).
"""
import pytest

import symox as sx
from symox import (S, Interval, Intersection, Complement, FiniteSet,
                   ProductSet, Range, pi, oo, sin, cos, tan, exp, log, sqrt,
                   Abs, I, Lambda, ImageSet, Rational, Symbol, srepr,
                   solveset, diophantine)
from symox.solvers.solveset import invert_complex, invert_real

x = Symbol('x')
r = Symbol('r', real=True)
p = Symbol('p', positive=True)
n = Symbol('n', integer=True)
m = Symbol('m', integer=True)
k = Symbol('k', integer=True, positive=True)

INVERT = [
    ('invert_complex(n*I, 1, n)',
     "(Symbol('n', integer=True), FiniteSet(Mul(Integer(-1), I)))"),
    ('invert_complex(2*n, 1, n)',
     "(Symbol('n', integer=True), FiniteSet(Rational(1, 2)))"),
    ('invert_real(2*n, 1, n)',
     "(Symbol('n', integer=True), FiniteSet(Rational(1, 2)))"),
    ('invert_real(n**2, 2, n)',
     "(Symbol('n', integer=True), FiniteSet(Pow(Integer(2), Rational(1, 2)), Mul(Integer(-1), Pow(Integer(2), Rational(1, 2)))))"),
    ('invert_real(p, -1, p)',
     "(Symbol('p', positive=True), FiniteSet(Integer(-1)))"),
    ('invert_real(2*p, -1, p)',
     "(Symbol('p', positive=True), FiniteSet(Rational(-1, 2)))"),
    ('invert_complex(p*I, 1, p)',
     "(Symbol('p', positive=True), FiniteSet(Mul(Integer(-1), I)))"),
    ('invert_real(k, Rational(1,2), k)',
     "(Symbol('k', integer=True, positive=True), FiniteSet(Rational(1, 2)))"),
    ('invert_real(k+1, 0, k)',
     "(Symbol('k', integer=True, positive=True), FiniteSet(Integer(-1)))"),
    ('invert_complex(2*n*I, 3, n)',
     "(Symbol('n', integer=True), FiniteSet(Mul(Integer(-1), Rational(3, 2), I)))"),
    ('invert_complex(r*I, 1, r)',
     "(Symbol('r', real=True), FiniteSet(Mul(Integer(-1), I)))"),
    ('invert_complex(r, I, r)',
     "(Symbol('r', real=True), FiniteSet(I))"),
    ('invert_real(p**2, 4, p)',
     "(Symbol('p', positive=True), FiniteSet(Integer(-2), Integer(2)))"),
    ('invert_real(k**2, 4, k)',
     "(Symbol('k', integer=True, positive=True), FiniteSet(Integer(-2), Integer(2)))"),
    ('invert_complex(p+I, 1, p)',
     "(Symbol('p', positive=True), FiniteSet(Add(Integer(1), Mul(Integer(-1), I))))"),
    ('invert_real(2*k, 3, k)',
     "(Symbol('k', integer=True, positive=True), FiniteSet(Rational(3, 2)))"),
]

DIOPHANTINE = [
    ('diophantine(n**2 - m, syms=(n, m), permute=True)',
     "['(t, t**2)']"),
    ('diophantine(n**2 + n - m, syms=(n, m), permute=True)',
     "['(t, t**2 + t)']"),
    ('diophantine(2*m - n**2, syms=(n, m), permute=True)',
     "['(-2*t, 2*t**2)']"),
    ('diophantine(n**2 - 2*m, syms=(n, m), permute=True)',
     "['(-2*t, 2*t**2)']"),
    ('diophantine(n**2 + n - 2*m, syms=(n, m), permute=True)',
     "['(2*t + 1, 2*t**2 + 3*t + 1)', '(2*t, 2*t**2 + t)']"),
    ('diophantine(4*n**2 + 4*n*m + m**2 + n + 1, syms=(n, m))',
     "['(-t**2 - 1, 2*t**2 - t + 2)']"),
    ('diophantine(n**2 - 2*n*m + m**2 - n, syms=(n, m))',
     "['(t**2, t**2 + t)']"),
    ('diophantine(n**2 + 2*n*m + m**2 - 3*n - 3*m, syms=(n, m))',
     "['(-t_0, t_0)', '(3 - t_0, t_0)']"),
    ('diophantine(2*pi*n - pi*m, syms=(n, m), permute=True)',
     "['(t_0, 2*t_0)']"),
    ('diophantine(pi*n - pi*m, syms=(n, m), permute=True)',
     "['(t_0, t_0)']"),
    ('diophantine(n**2 + m, syms=(n, m))',
     "['(t, -t**2)']"),
    ('diophantine(n**2 - 4*m - 1, syms=(n, m))',
     "['(1 - 4*t, 4*t**2 - 2*t)', '(3 - 4*t, 4*t**2 - 6*t + 2)']"),
    ('diophantine((n - m)**2 - n - m, syms=(n, m))',
     "['(2*t**2 + t, 2*t**2 - t)', '(2*t**2 - t, 2*t**2 - 3*t + 1)']"),
    ('diophantine(m - n**2, syms=(m, n))',
     "['(t**2, t)']"),
    ('diophantine(n**2 - m - 3, syms=(n, m), permute=True)',
     "['(t, t**2 - 3)']"),
    ('diophantine(3*n**2 - m, syms=(n, m))',
     "['(t, 3*t**2)']"),
    ('diophantine(n**2 + 2*n + 1 - m, syms=(n, m))',
     "['(t, t**2 + 2*t + 1)']"),
    ('diophantine(n**2 + 2*n*m + m**2 - 1, syms=(n, m))',
     "['(-t_0 - 1, t_0)', '(1 - t_0, t_0)']"),
    ('diophantine(x**2 - 2*x*m + m**2)',
     "['(t_0, t_0)']"),
]

SETS = [
    ('Intersection(S.Integers, ImageSet(Lambda(n, n**2), S.Integers))',
     "ImageSet(Lambda(Tuple(Symbol('n', integer=True)), Pow(Symbol('n', integer=True), Integer(2))), Integers)"),
    ('Intersection(S.Integers, ImageSet(Lambda(n, n**2 + n), S.Integers))',
     "ImageSet(Lambda(Tuple(Symbol('n', integer=True)), Add(Pow(Symbol('n', integer=True), Integer(2)), Symbol('n', integer=True))), Integers)"),
    ('Intersection(ImageSet(Lambda(n, n**2), S.Integers), ImageSet(Lambda(n, 2*n), S.Integers))',
     "ImageSet(Lambda(Tuple(Symbol('n', integer=True)), Mul(Integer(4), Pow(Symbol('n', integer=True), Integer(2)))), Integers)"),
    ('Intersection(ImageSet(Lambda(n, 2*n*pi), S.Integers), ImageSet(Lambda(n, n*pi), S.Integers))',
     "ImageSet(Lambda(Tuple(Symbol('n', integer=True)), Mul(Integer(2), pi, Symbol('n', integer=True))), Integers)"),
    ('Intersection(ImageSet(Lambda(n, n*pi), S.Integers), ImageSet(Lambda(n, 2*n*pi + pi), S.Integers))',
     "ImageSet(Lambda(Tuple(Symbol('n', integer=True)), Add(Mul(Integer(2), pi, Symbol('n', integer=True)), pi)), Integers)"),
    ('Intersection(S.Integers, ImageSet(Lambda(n, 2*n**2 + 1), S.Integers))',
     "ImageSet(Lambda(Tuple(Symbol('n', integer=True)), Add(Mul(Integer(2), Pow(Symbol('n', integer=True), Integer(2))), Integer(1))), Integers)"),
    ('Intersection(ImageSet(Lambda(n, 3*n), S.Integers), ImageSet(Lambda(n, n**2), S.Integers))',
     "ImageSet(Lambda(Tuple(Symbol('n', integer=True)), Mul(Integer(9), Pow(Symbol('n', integer=True), Integer(2)))), Integers)"),
    ('Range(0, 10, 2).contains(2*n)',
     "And(GreaterThan(Mul(Integer(2), Symbol('n', integer=True)), Integer(0)), LessThan(Mul(Integer(2), Symbol('n', integer=True)), Integer(8)))"),
    ('Range(0, 10, 2).contains(2*n + 1)',
     'false'),
    ('Range(1, 10, 3).contains(3*n + 1)',
     "And(GreaterThan(Add(Mul(Integer(3), Symbol('n', integer=True)), Integer(1)), Integer(1)), LessThan(Add(Mul(Integer(3), Symbol('n', integer=True)), Integer(1)), Integer(7)))"),
    ('Range(0, oo, 2).contains(2*n)',
     "GreaterThan(Mul(Integer(2), Symbol('n', integer=True)), Integer(0))"),
    ('Range(-oo, 0, 2).contains(2*n)',
     "LessThan(Mul(Integer(2), Symbol('n', integer=True)), Integer(-2))"),
    ('Range(10, 0, -2).contains(2*n)',
     "And(GreaterThan(Mul(Integer(2), Symbol('n', integer=True)), Integer(2)), LessThan(Mul(Integer(2), Symbol('n', integer=True)), Integer(10)))"),
    ('Range(0, 10, 2).contains(4*n)',
     "And(GreaterThan(Mul(Integer(4), Symbol('n', integer=True)), Integer(0)), LessThan(Mul(Integer(4), Symbol('n', integer=True)), Integer(8)))"),
    ('Range(0, 10, 2).contains(2*n**2)',
     "LessThan(Mul(Integer(2), Pow(Symbol('n', integer=True), Integer(2))), Integer(8))"),
    ('Range(0, 3, 5).contains(n)',
     "Equality(Symbol('n', integer=True), Integer(0))"),
    ('Range(0, 3, 5).contains(5*n)',
     "Equality(Mul(Integer(5), Symbol('n', integer=True)), Integer(0))"),
    ('Range(0, 10, 2).contains(2*n + 2*m)',
     "And(GreaterThan(Add(Mul(Integer(2), Symbol('m', integer=True)), Mul(Integer(2), Symbol('n', integer=True))), Integer(0)), LessThan(Add(Mul(Integer(2), Symbol('m', integer=True)), Mul(Integer(2), Symbol('n', integer=True))), Integer(8)))"),
    ('Range(1, 10, 2).contains(2*n + 1)',
     "And(GreaterThan(Add(Mul(Integer(2), Symbol('n', integer=True)), Integer(1)), Integer(1)), LessThan(Add(Mul(Integer(2), Symbol('n', integer=True)), Integer(1)), Integer(9)))"),
    ('ImageSet(Lambda(n, 2*n), S.Integers).contains(x)',
     "Contains(Symbol('x'), ImageSet(Lambda(Tuple(Symbol('n', integer=True)), Mul(Integer(2), Symbol('n', integer=True))), Integers))"),
    ('ImageSet(Lambda(n, 2*n + 1), S.Integers).contains(x)',
     "Contains(Symbol('x'), ImageSet(Lambda(Tuple(Symbol('n', integer=True)), Add(Mul(Integer(2), Symbol('n', integer=True)), Integer(1))), Integers))"),
    ('ImageSet(Lambda(n, 2*n), S.Integers).contains(m)',
     "Contains(Symbol('m', integer=True), ImageSet(Lambda(Tuple(Symbol('n', integer=True)), Mul(Integer(2), Symbol('n', integer=True))), Integers))"),
    ('ImageSet(Lambda(n, n*pi), S.Integers).contains(x)',
     "Contains(Symbol('x'), ImageSet(Lambda(Tuple(Symbol('n', integer=True)), Mul(pi, Symbol('n', integer=True))), Integers))"),
    ('ProductSet(Interval(0, 1), Interval(0, 1)).contains((x, r))',
     "And(LessThan(Integer(0), Symbol('r', real=True)), LessThan(Integer(0), Symbol('x')), LessThan(Symbol('r', real=True), Integer(1)), LessThan(Symbol('x'), Integer(1)))"),
    ('ProductSet(Interval(0, 1), Interval(0, 1)).contains((r, 2))',
     'false'),
    ('ProductSet(Interval(0, 1), Interval(0, 1)).contains((Rational(1,2), 2))',
     'false'),
    ('ProductSet(Interval(0, 1), Interval(0, 1)).contains((Rational(1,2), Rational(1,3)))',
     'true'),
    ('ProductSet(Interval(0, 1), Interval(0, 1)).contains((x, r, 1))',
     'false'),
    ('ProductSet(S.Integers, S.Reals).contains((n, r))',
     'true'),
    ('ProductSet(S.Integers, S.Reals).contains((x, r))',
     "Contains(Symbol('x'), Integers)"),
    ('ProductSet(FiniteSet(1, 2), Interval(0, 1)).contains((x, 0))',
     "Or(Equality(Integer(1), Symbol('x')), Equality(Integer(2), Symbol('x')))"),
    ('solveset(n*I - 1, n)',
     'FiniteSet(Mul(Integer(-1), I))'),
    ('solveset(2*n - 1, n)',
     'FiniteSet(Rational(1, 2))'),
    ('solveset(n**2 - 2, n, S.Reals)',
     'FiniteSet(Pow(Integer(2), Rational(1, 2)), Mul(Integer(-1), Pow(Integer(2), Rational(1, 2))))'),
]


@pytest.mark.parametrize('call,expected', INVERT)
def test_invert_ignores_symbol_assumptions(call, expected):
    got = eval(call)
    assert '(' + ', '.join(srepr(i) for i in got) + ')' == expected


@pytest.mark.parametrize('call,expected', DIOPHANTINE)
def test_diophantine_reference_families(call, expected):
    assert repr(sorted(str(t) for t in eval(call))) == expected


@pytest.mark.parametrize('call,expected', SETS)
def test_sets_membership_and_intersections(call, expected):
    assert srepr(eval(call)) == expected


def test_diophantine_rejects_irrational_coefficients():
    with pytest.raises(TypeError):
        diophantine(sqrt(2)*n - m, syms=(n, m))


# The reference misses solutions here; symox keeps the complete sets.
CORRECTED = [
    # SymPy 1.14 gives {0}: exp(I*x) = 1 for every real x = 2*n*pi
    ('solveset(exp(I*x) - 1, x, S.Reals)', 'ImageSet'),
    # SymPy 1.14 gives {pi/2}: 5*pi/2 < 10 is also a solution
    ('solveset(exp(I*x) - I, x, Interval(0, 10))',
     'FiniteSet(Mul(Rational(1, 2), pi), Mul(Rational(5, 2), pi))'),
]


@pytest.mark.parametrize('call,expected', CORRECTED)
def test_complete_periodic_solutions_kept(call, expected):
    assert srepr(eval(call)).startswith(expected)
