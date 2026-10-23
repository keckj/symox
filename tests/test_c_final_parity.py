"""Regression tests for the last undocumented "c" gaps of the SymPy 1.14
doc-example audit (harmonic, loggamma, solve, lambdify, autowrap).

Expected values are hard-coded from SymPy 1.14."""
import warnings

import pytest

from symox import (
    Abs, Eq, Function, Piecewise, PoleError, Rational, S, cos, exp, expand_func,
    harmonic, lambdify, loggamma, reduce_inequalities, series, simplify, sin, solve,
    sqrt, symbols,
)

x, y, n = symbols('x y n')


# --- harmonic at rational arguments ------------------------------------------

@pytest.mark.parametrize('arg, expected', [
    (S(1)/3, '-log(6) - sqrt(3)*pi/6 + 2*Sum(log(sin(_k*pi/3))*cos(2*_k*pi/3), (_k, 1, 1))'
             ' + 3*Sum(1/(3*_k + 1), (_k, 0, 0))'),
    (S(25)/7, '-log(14) + pi*cot(3*pi/7)/2 + 2*Sum(log(sin(_k*pi/7))*cos(8*_k*pi/7), (_k, 1, 3))'
              ' + 7*Sum(1/(7*_k + 4), (_k, 0, 3))'),
    (S(1)/2, '-log(4) + 2*Sum(log(sin(_k*pi/2))*cos(_k*pi), (_k, 1, 0))'
             ' + 2*Sum(1/(2*_k + 1), (_k, 0, 0))'),
])
def test_harmonic_rational_expand_func(arg, expected):
    assert str(harmonic(arg)) == 'harmonic(%s)' % arg
    assert str(expand_func(harmonic(arg))) == expected


def test_harmonic_rational_doit_values():
    assert str(expand_func(harmonic(S(1)/3)).doit()) == '-log(6) - sqrt(3)*pi/6 - log(sqrt(3)/2) + 3'
    assert str(expand_func(harmonic(S(1)/2)).doit()) == '2 - log(4)'
    assert str(expand_func(harmonic(S(7)/2)).doit()) == '352/105 - log(4)'
    he = expand_func(harmonic(S(25)/7)).doit()
    assert str(he.n(40)) == '1.983697455232980674869851942390639915940'
    assert str(harmonic(S(25)/7).n(40)) == '1.983697455232980674869851942390639915940'
    assert str(simplify(he).n(40)) == '1.983697455232980674869851942390639915940'


def test_harmonic_expand_func_left_alone():
    # negative rationals and orders m != 1 are left as they are
    assert str(expand_func(harmonic(S(-1)/3))) == 'harmonic(-1/3)'
    assert str(expand_func(harmonic(S(1)/3, 2))) == 'harmonic(1/3, 2)'


def test_harmonic_integer_offsets():
    assert str(expand_func(harmonic(n - 4))) == 'harmonic(n) - 1/(n - 1) - 1/(n - 2) - 1/(n - 3) - 1/n'
    assert str(expand_func(harmonic(n + 4))) == 'harmonic(n) + 1/(n + 4) + 1/(n + 3) + 1/(n + 2) + 1/(n + 1)'
    assert str(expand_func(harmonic(n - 1))) == 'harmonic(n) - 1/n'


# --- loggamma ----------------------------------------------------------------

@pytest.mark.parametrize('arg, expected', [
    (S(16)/3, '-5*log(3) + loggamma(1/3) + Sum(log(3*_k - 2), (_k, 1, 5))'),
    (S(-5)/3, '-Sum(log(3*_k - 1), (_k, 1, 2)) + loggamma(1/3) + 2*log(3) - 2*I*pi'),
    (S(-1)/3, '-Sum(log(3*_k - 2), (_k, 1, 1)) + loggamma(2/3) + log(3) - I*pi'),
    (S(2)/3, 'loggamma(2/3)'),
])
def test_loggamma_rational_expand_func(arg, expected):
    assert str(expand_func(loggamma(arg))) == expected


def test_loggamma_rational_expand_func_values():
    assert str(expand_func(loggamma(S(16)/3)).doit()) == \
        '-5*log(3) + loggamma(1/3) + log(4) + log(7) + log(10) + log(13)'
    for value in (S(19)/4, S(-11)/4):
        assert abs(expand_func(loggamma(value)).n(30) - loggamma(value).n(30)) < 1e-25


@pytest.mark.parametrize('expr, order, expected', [
    (loggamma(x), 4, '-log(x) - EulerGamma*x + pi**2*x**2/12 - x**3*zeta(3)/3 + O(x**4)'),
    (loggamma(2*x), 3, '-log(2) - log(x) - 2*EulerGamma*x + pi**2*x**2/3 + O(x**3)'),
    (loggamma(x**2 + x), 3, '-log(x) + x*(-1 - EulerGamma) + x**2*(-EulerGamma + 1/2 + pi**2/12) + O(x**3)'),
    (loggamma(sin(x)), 3, '-log(x) - EulerGamma*x + x**2*(1/6 + pi**2/12) + O(x**3)'),
    (loggamma(x + 1), 4, '-EulerGamma*x + pi**2*x**2/12 - x**3*zeta(3)/3 + O(x**4)'),
])
def test_loggamma_series_at_zero(expr, order, expected):
    assert str(series(expr, x, 0, order)) == expected


def test_loggamma_series_doc_example_and_pole():
    assert str(series(loggamma(x), x, 0, 4).cancel()) == \
        '-log(x) - EulerGamma*x + pi**2*x**2/12 - x**3*zeta(3)/3 + O(x**4)'
    assert str(loggamma(x).nseries(x, 0, 3)) == '-log(x) - EulerGamma*x + pi**2*x**2/12 + O(x**3)'
    with pytest.raises(PoleError):
        series(loggamma(x - 1), x, 0, 3)


# --- solve: radical equations, implicit=True ---------------------------------

def test_solve_radical_equations_without_symbols():
    eq = x**Rational(1, 3) - x**Rational(1, 5) + Rational(1, 7)
    expected = '[CRootOf(7*x**5 - 7*x**3 + 1, 1)**15, CRootOf(7*x**5 - 7*x**3 + 1, 2)**15]'
    assert str(solve(eq)) == expected
    assert str(solve(eq, x)) == expected
    # the candidate 1/3 of unrad does not satisfy the equation
    assert solve((x**3 - 3*x**2)**Rational(1, 3) + 1 - x) == []
    assert solve(sqrt(x + 1) - sqrt(x) - 1) == [0]
    assert solve(sqrt(x + 1) - sqrt(x) - 1, dict=True) == [{x: 0}]


@pytest.mark.parametrize('args, kwargs, expected', [
    ((x + exp(x), x), {}, '[-exp(x)]'),
    ((x + exp(x),), {}, '[-exp(x)]'),
    ((x + exp(x)*y, x), {}, '[-y*exp(x)]'),
    ((x - y*exp(x) - 1, x), {}, '[y*exp(x) + 1]'),
    ((x**2 + y*x - x, x), {}, '[-x**2/(y - 1)]'),
    ((x + sin(x) - 1, x), {}, '[1 - sin(x)]'),
    ((x - y*exp(x) - 1, x), {'dict': True}, '[{x: y*exp(x) + 1}]'),
    ((x - y - 1, x), {}, '[y + 1]'),
    ((Eq(x, exp(-x)), x), {}, '[exp(-x)]'),
    (([x + exp(y), y - x*exp(x)], [x, y]), {}, '{x: -exp(y), y: -exp(x)*exp(y)}'),
    (([x - y - 1, Function('f')(x)*y], [x, y]), {}, '{x: 1, y: 0}'),
])
def test_solve_implicit(args, kwargs, expected):
    assert str(solve(*args, implicit=True, **kwargs)) == expected


def test_solve_implicit_false_is_plain_solve():
    assert solve(x - y - 1, x, implicit=False) == [y + 1]


# --- reduce_inequalities over several symbols --------------------------------

@pytest.mark.parametrize('ineqs, syms, expected', [
    ([0 <= x - 1, y < 3], [x, y], '(-oo < y) & (1 <= x) & (x < oo) & (y < 3)'),
    ([x**2 < 4, y > 1], [x, y], '(-2 < x) & (1 < y) & (x < 2) & (y < oo)'),
    ([x**2 < 4, y > 1], [], '(-2 < x) & (1 < y) & (x < 2) & (y < oo)'),
    ([x > 1, y < 3], [x], '(-oo < y) & (1 < x) & (x < oo) & (y < 3)'),
    ([Abs(x) < 4, y > 1], [x, y], '(-4 < x) & (1 < y) & (x < 4) & (y < oo)'),
    ([x > 1, x < 5, y**2 >= 4], [x, y],
     '(1 < x) & (x < 5) & (((2 <= y) & (y < oo)) | ((-oo < y) & (y <= -2)))'),
])
def test_reduce_inequalities_several_symbols(ineqs, syms, expected):
    assert str(reduce_inequalities(ineqs, syms)) == expected


def test_reduce_inequalities_coupled_symbols_raise():
    with pytest.raises(NotImplementedError) as info:
        reduce_inequalities([x > 1, y < 3, x + y > 0], [x, y])
    assert str(info.value) == '\ninequality has more than one symbol of interest.'


# --- lambdify namespace ------------------------------------------------------

def test_lambdify_globals_numpy():
    numpy = pytest.importorskip('numpy')
    f = lambdify(x, sin(x) + cos(x))
    assert f.__globals__['sin'] is numpy.sin
    assert f.__globals__['cos'] is numpy.cos
    assert f.__globals__['x'] == x
    assert f.__globals__['range'] is range
    assert abs(f(1.0) - 1.3817732906760363) < 1e-15


def test_lambdify_globals_modules():
    import math
    numpy = pytest.importorskip('numpy')
    assert lambdify(x, sin(x), 'math').__globals__['sin'] is math.sin
    assert lambdify(x, sin(x), 'math').__globals__['ceiling'] is math.ceil

    def mysin(value):
        return numpy.sin(value)
    k = lambdify(x, sin(x), [{'sin': mysin}, 'numpy'])
    assert k.__globals__['sin'] is mysin
    assert k.__globals__['cos'] is numpy.cos
    assert abs(k(1) - 0.8414709848078965) < 1e-15


def test_testing_pytest_ignore_warnings():
    numpy = pytest.importorskip('numpy')
    from symox.testing.pytest import ignore_warnings
    f = lambdify(x, Piecewise((x, x <= 1), (1/x, x > 1)), 'numpy')
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        with ignore_warnings(RuntimeWarning):
            values = f(numpy.array([-1, 0, 1, 2]))
    assert list(values) == [-1.0, 0.0, 1.0, 0.5]


# --- autowrap helpers --------------------------------------------------------

def test_autowrap_single_helper_tuple_is_normalized(monkeypatch):
    from symox.utilities import autowrap as module
    seen = []

    class Recorder:
        def __init__(self, *args, **kwargs):
            pass

        def wrap_code(self, routine, helpers):
            seen.append([h.name for h in helpers])
            return routine

    t = symbols('t')
    h = Function('helper_func')
    monkeypatch.setitem(module.autowrap.__globals__, 'DummyWrapper', Recorder)
    module.autowrap(3*x + h(t), args=[x, t], backend='dummy', helpers=('helper_func', 4*x, [x]))
    module.autowrap(3*x + h(t), args=[x, t], backend='dummy', helpers=[('helper_func', 4*x, [x])])
    module.autowrap(3*x, args=[x], backend='dummy')
    assert seen == [['helper_func'], ['helper_func'], []]
