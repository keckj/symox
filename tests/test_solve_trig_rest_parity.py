"""Trig equations in solve and solveset, values from SymPy 1.14."""
import pytest
from symox import (I, S, Symbol, symbols, solve, solveset, sin, cos, tan, cot,
                   sec, csc, sinh, exp, pi, sqrt)

x, y, a = symbols('x y a')
xr = Symbol('xr', real=True)
p = Symbol('p', positive=True)

CASES = [
    ('solve(sin(x)+cos(x), x)',
     '[-pi/4]'),
    ('solve(sin(x)+cos(x)-1, x)',
     '[0, pi/2]'),
    ('solve(sin(2*x)-cos(x), x)',
     '[-pi/2, pi/6, pi/2, 5*pi/6]'),
    ('solve(sin(x)**2-cos(x)**2, x)',
     '[-pi/4, pi/4]'),
    ('solve(sin(x)-cos(x), x)',
     '[pi/4]'),
    ('solve(cos(2*x)-cos(x), x)',
     '[0, -2*pi/3, 2*pi/3]'),
    ('solve(sin(x)*cos(x)-S(1)/4, x)',
     '[2*atan(-2 + 2*sqrt(2 - sqrt(3)) + sqrt(3)), -2*atan(-sqrt(3) + 2*sqrt(2 - sqrt(3)) + 2), -2*atan(sqrt(3) + 2 + 2*sqrt(sqrt(3) + 2)), -2*atan(-2*sqrt(sqrt(3) + 2) + sqrt(3) + 2)]'),
    ('solve(tan(x)-cot(x), x)',
     '[-pi/4, pi/4]'),
    ('solve(sin(x)/cos(x)-1, x)',
     '[pi/4]'),
    ('solve(3*sin(x)+4*cos(x), x)',
     '[-atan(4/3)]'),
    ('solve(3*sin(x)+4*cos(x)-5, x)',
     '[2*atan(1/3)]'),
    ('solve(sin(x)+cos(x), x, dict=True)',
     '[{x: -pi/4}]'),
    ('solve(sin(x)+cos(x)-1, x, dict=True)',
     '[{x: 0}, {x: pi/2}]'),
    ('solve(sin(x)+cos(x), x, check=False)',
     '[-pi/4]'),
    ('solve(sin(x)+cos(x)-a, x)',
     '[-2*atan((sqrt(2 - a**2) - 1)/(a + 1)), 2*atan((sqrt(2 - a**2) + 1)/(a + 1))]'),
    ('solve(a*sin(x)+cos(x), x)',
     '[-atan(1/a)]'),
    ('solve(sin(xr)+cos(xr), xr)',
     '[-pi/4]'),
    ('solve(sin(p)+cos(p), p)',
     '[]'),
    ('solve(sin(x)**2+cos(x)-1, x)',
     '[0, -pi/2, pi/2]'),
    ('solve(sin(x)+sin(2*x), x)',
     '[0, -2*pi/3, 2*pi/3]'),
    ('solve(sin(x)+cos(x)+1, x)',
     '[-pi/2]'),
    ('solve(sin(x)-cos(2*x), x)',
     '[-pi/2, pi/6, 5*pi/6]'),
    ('solve(sin(x)+I*cos(x), x)',
     '[]'),
    ('solve([sin(x)+cos(x), y-1], [x, y])',
     '[(-pi/4, 1)]'),
    ('solve(sin(x)**2+2*sin(x)*cos(x)-3*cos(x)**2, x)',
     '[-3*pi/4, pi/4, 2*atan(1/3 - sqrt(10)/3), 2*atan(1/3 + sqrt(10)/3)]'),
    ('solve(sqrt(3)*sin(x)+cos(x)-1, x)',
     '[0, 2*pi/3]'),
    ('solveset(sin(a*x)-1, x)',
     'ConditionSet(x, Ne(a, 0), ImageSet(Lambda(_n, (2*_n*pi + pi/2)/a), Integers))'),
    ('solveset(cos(a*x), x)',
     'ConditionSet(x, Ne(a, 0), Union(ImageSet(Lambda(_n, (2*_n*pi + pi/2)/a), Integers), ImageSet(Lambda(_n, (2*_n*pi + 3*pi/2)/a), Integers)))'),
    ('solveset(sin(a*x)-1, x, S.Reals)',
     'ConditionSet(x, Ne(a, 0), ImageSet(Lambda(_n, (2*_n*pi + pi/2)/a), Integers))'),
    ('solveset(sin(2*a*x), x)',
     'ConditionSet(x, Ne(2*a, 0), Union(ImageSet(Lambda(_n, (2*_n*pi + pi)/(2*a)), Integers), ImageSet(Lambda(_n, _n*pi/a), Integers)))'),
    ('solveset(tan(a*x)-1, x)',
     'ConditionSet(x, Ne(a, 0), ImageSet(Lambda(_n, (_n*pi + pi/4)/a), Integers))'),
    ('solveset(sinh(a*x)-1, x)',
     'ConditionSet(x, Ne(a, 0), Union(ImageSet(Lambda(_n, (2*_n*I*pi + log(1 + sqrt(2)))/a), Integers), ImageSet(Lambda(_n, (2*_n*I*pi - log(1 + sqrt(2)) + I*pi)/a), Integers)))'),
    ('solve(tan(x)-sin(x), x)',
     '[0, -pi, pi, 2*pi]'),
    ('solve(sin(x)+cos(x)-sqrt(2), x)',
     '[pi/4]'),
    ('solve(sec(x)-2, x)',
     '[-pi/3, pi/3]'),
    ('solve(csc(x)-2, x)',
     '[pi/6, 5*pi/6]'),
    ('solve(sec(x)+1, x)',
     '[pi]'),
    ('solve(csc(x)+1, x)',
     '[-pi/2]'),
    ('solve(sec(x)-sqrt(2), x)',
     '[-pi/4, pi/4]'),
    ('solve(sec(x)-1, x)',
     '[0]'),
]


@pytest.mark.parametrize("call, expected", CASES)
def test_trig_parity(call, expected):
    assert str(eval(call)) == expected


def test_undecided_stays_unsolved():
    for e in (exp(x) - sin(x), x - sin(x)):
        with pytest.raises(NotImplementedError):
            solve(e, x)
    assert str(solveset(sin(a*x)**2 - 1, x)) ==         'ConditionSet(x, Eq(sin(a*x)**2 - 1, 0), Complexes)'
