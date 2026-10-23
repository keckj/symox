"""Regression tests for the second batch of solve/solveset parity fixes.

Expected strings are the outputs of the reference SymPy 1.14.
"""

import pytest

import symox
from symox import Symbol, pi, simplify, sin, solve, symbols

x, y, z, a, b, n = symbols('x y z a b n')
ap = Symbol('ap', positive=True)

CASES = [
    ('solve(sec(x) - 3)',
     '[I*(log(3) - log(1 - 2*sqrt(2)*I)), I*(log(3) - log(1 + 2*sqrt(2)*I))]'),
    ('solve(csc(x) - 5)',
     '[I*(log(5) - log(-2*sqrt(6) + I)), I*(log(5) - log(2*sqrt(6) + I))]'),
    ('solve(csc(x) - 3)',
     '[I*(log(3) - log(-2*sqrt(2) + I)), I*(log(3) - log(2*sqrt(2) + I))]'),
    ('solve(sec(x) + 3)',
     '[I*(log(3) - log(-1 - 2*sqrt(2)*I)), I*(log(3) - log(-1 + 2*sqrt(2)*I))]'),
    ('solve(sec(x) - Rational(1,2))',
     '[-I*log(2 - sqrt(3)), -I*log(sqrt(3) + 2)]'),
    ('solve(sec(x) - 2)',
     '[-pi/3, pi/3]'),
    ('solve(csc(x) - 2)',
     '[pi/6, 5*pi/6]'),
    ('solve(sec(x) + 1)',
     '[pi]'),
    ('solve(csc(x) - sqrt(2))',
     '[pi/4, 3*pi/4]'),
    ('solveset(exp(exp(x)) - 1, x)',
     'ImageSet(Lambda(((_k1, _n),), I*(2*_k1*pi + arg(_n*I)) + log(2*pi*Abs(_n))), ProductSet(Integers, Integers))'),
    ('solveset(exp(exp(x)) - 2, x)',
     'ImageSet(Lambda(((_k1, _n),), I*(2*_k1*pi + arg(2*_n*I*pi + log(2))) + log(Abs(2*_n*I*pi + log(2)))), ProductSet(Integers, Integers))'),
    ('solveset(exp(exp(2*x)) - 1, x)',
     'ImageSet(Lambda(x, x/2), ImageSet(Lambda(((_k1, _n),), I*(2*_k1*pi + arg(_n*I)) + log(2*pi*Abs(_n))), ProductSet(Integers, Integers)))'),
    ('solveset(exp(exp(x) + 1) - 1, x)',
     'ImageSet(Lambda(((_k1, _n),), I*(2*_k1*pi + arg(2*_n*I*pi - 1)) + log(Abs(2*_n*I*pi - 1))), ProductSet(Integers, Integers))'),
    ('solveset(exp(exp(x)) - 1, x, S.Reals)',
     'EmptySet'),
    ('solveset(exp(exp(exp(x))) - 1, x)',
     'ImageSet(Lambda(((_k2, _k1, _n),), I*(2*_k2*pi + arg(I*(2*_k1*pi + arg(_n*I)) + log(2*pi*Abs(_n)))) + log(Abs(I*(2*_k1*pi + arg(_n*I)) + log(2*pi*Abs(_n))))), ProductSet(Integers, Integers, Integers))'),
    ('solveset(exp(2*exp(x)) - 3, x)',
     'ImageSet(Lambda(((_k1, _n),), I*(2*_k1*pi + arg(_n*I*pi + log(3)/2)) + log(Abs(_n*I*pi + log(3)/2))), ProductSet(Integers, Integers))'),
    ('solveset(exp(exp(x)) - E, x)',
     'ImageSet(Lambda(((_k1, _n),), I*(2*_k1*pi + arg(2*_n*I*pi + 1)) + log(Abs(2*_n*I*pi + 1))), ProductSet(Integers, Integers))'),
    ('solveset(cos(a*x) - b, x)',
     'ConditionSet(x, Ne(a, 0), Union(ImageSet(Lambda(_n, (2*_n*pi - acos(b))/a), Integers), ImageSet(Lambda(_n, (2*_n*pi + acos(b))/a), Integers)))'),
    ('solveset(sin(a*x) - b, x)',
     'ConditionSet(x, Ne(a, 0), Union(ImageSet(Lambda(_n, (2*_n*pi + asin(b))/a), Integers), ImageSet(Lambda(_n, (2*_n*pi - asin(b) + pi)/a), Integers)))'),
    ('solveset(cos(a*x) - b, x, S.Reals)',
     'ConditionSet(x, (-1 <= b) & (b <= 1) & Ne(a, 0), Union(ImageSet(Lambda(_n, (2*_n*pi - acos(b))/a), Integers), ImageSet(Lambda(_n, (2*_n*pi + acos(b))/a), Integers)))'),
    ('solveset(sin(a*x) - cos(a*x), x)',
     'ConditionSet(x, Ne(a, 0), Union(ImageSet(Lambda(_n, (2*_n*pi - 3*pi/4)/a), Integers), ImageSet(Lambda(_n, (2*_n*pi + pi/4)/a), Integers)))'),
    ('solveset(tan(a*x) - b, x)',
     'ConditionSet(x, Ne(a, 0), ImageSet(Lambda(_n, (_n*pi + atan(b))/a), Integers))'),
    ('solveset(sin(a*x)*cos(a*x), x)',
     'Union(ConditionSet(x, Ne(a, 0), Union(ImageSet(Lambda(_n, (2*_n*pi + pi/2)/a), Integers), ImageSet(Lambda(_n, (2*_n*pi + 3*pi/2)/a), Integers))), ConditionSet(x, Ne(a, 0), Union(ImageSet(Lambda(_n, (2*_n*pi + pi)/a), Integers), ImageSet(Lambda(_n, 2*_n*pi/a), Integers))))'),
    ('solveset(cos(a*x) - Rational(1,2), x)',
     'ConditionSet(x, Ne(a, 0), Union(ImageSet(Lambda(_n, (2*_n*pi + pi/3)/a), Integers), ImageSet(Lambda(_n, (2*_n*pi + 5*pi/3)/a), Integers)))'),
    ('solveset(sin(a*x + 1) - b, x)',
     'ConditionSet(x, Ne(a, 0), Union(ImageSet(Lambda(_n, (2*_n*pi + asin(b) - 1)/a), Integers), ImageSet(Lambda(_n, (2*_n*pi - asin(b) - 1 + pi)/a), Integers)))'),
    ('nonlinsolve([x**2 + 1/y - 2, x - y], [x, y])',
     '{(1, 1), (-1/2 + sqrt(5)/2, -1/2 + sqrt(5)/2), (-sqrt(5)/2 - 1/2, -sqrt(5)/2 - 1/2)}'),
    ('solve(sec(x)**2 - 9)',
     '[-I*log(-1/3 - 2*sqrt(2)*I/3), -I*log(-1/3 + 2*sqrt(2)*I/3), -I*log(1/3 - 2*sqrt(2)*I/3), -I*log(1/3 + 2*sqrt(2)*I/3)]'),
]


@pytest.mark.parametrize("call, expected", CASES)
def test_reference_output(call, expected):
    namespace = dict(vars(symox), x=x, y=y, a=a, b=b, ap=ap)
    assert str(eval(call, namespace)) == expected


def test_sin3x_keeps_pi():
    # Accepted divergence: symox also returns pi, which is a genuine root.
    sols = solve(sin(3*x) - sin(x))
    assert pi in sols
    for s in sols:
        assert simplify((sin(3*x) - sin(x)).subs(x, s)) == 0
