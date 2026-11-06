"""The accumulation bounds arithmetic and the oscillating limit are native
(``symox::limits::accum``): the Python modules only admit the operands and
route to ``_core``. Values are those of SymPy 1.14."""
import pytest

import symox
from symox import (E, FiniteSet, I, Limit, Rational, Symbol, cos, exp, frac,
                   limit, log, oo, pi, sin, tan, zoo)
from symox import _core
from symox.calculus import _accum_eval, _accum_limit
from symox.calculus import accumulationbounds as abmod
from symox.calculus.accumulationbounds import AccumBounds as AB

x = Symbol('x')


def test_python_modules_hold_no_bounds_algorithm():
    for name in ('_walk', '_combine', '_bounded_value', '_native', '_bounded_affine',
                 '_held_oscillation', '_real_infinity_times_sign', '_function'):
        assert not hasattr(_accum_limit, name)
    for name in ('_sin', '_cos', '_tan', '_cot', '_exp', '_log', '_frac', '_hyperbolic',
                 '_EVAL', '_real_infinity'):
        assert not hasattr(_accum_eval, name)
    for name in ('_flat_mul', '_extended_real'):
        assert not hasattr(abmod, name)
    for name in ('_oscillating_limit', '_accum_evaluate', '_accum_op', '_accum_new',
                 '_accum_fold', '_accum_is_ge', '_accum_members', '_accum_intersection',
                 '_accum_union', '_accum_contains', '_accum_unary'):
        assert callable(getattr(_core, name))


class _Spy:
    def __init__(self, monkeypatch, name):
        self.calls = 0
        original = getattr(_core, name)

        def spy(*args):
            self.calls += 1
            return original(*args)
        monkeypatch.setattr(_core, name, spy)


@pytest.mark.parametrize('e, z0, expected', [
    (sin(x), oo, AB(-1, 1)),
    (sin(1/x), 0, AB(-1, 1)),
    (exp(sin(x)), oo, AB(exp(-1), E)),
    (sin(x)**2, oo, AB(0, 1)),
    (1/(2 + sin(x)), oo, AB(Rational(1, 3), 1)),
    (frac(x), oo, AB(0, 1)),
    (tan(x), oo, AB(-oo, oo)),
    (sin(x) + cos(x), oo, AB(-2, 2)),
    (frac(x)**(1/x), oo, AB(0, 1)),
    (log(2 + sin(x)), oo, AB(0, log(3))),
])
def test_oscillating_limit_is_native(monkeypatch, e, z0, expected):
    spy = _Spy(monkeypatch, '_oscillating_limit')
    assert limit(e, x, z0) == expected
    assert spy.calls >= 1


def test_two_sided_and_held():
    assert limit(sin(1/x), x, 0, '+-') == AB(-1, 1)
    # the reference keeps an oscillating base under an infinite exponent
    assert isinstance(limit(sin(x)**x, x, oo), Limit)
    assert limit(sin(x)/x, x, oo) == 0
    assert limit(cos(x) + x, x, oo) == oo
    # log of bounds reaching zero or below is nan: the reference holds
    assert isinstance(limit(log(sin(x)), x, oo), Limit)
    assert isinstance(limit(log(sin(x) + Rational(1, 2)), x, oo), Limit)
    assert limit(log(1 + sin(x)), x, oo) == AB(-oo, log(2))


def test_arithmetic_is_native(monkeypatch):
    spy = _Spy(monkeypatch, '_accum_op')
    assert AB(1, 2) + 1 == AB(2, 3)
    assert AB(1, 2)*-2 == AB(-4, -2)
    assert 1/AB(1, 2) == AB(Rational(1, 2), 1)
    assert AB(-1, 1)**2 == AB(0, 1)
    assert AB(1, 2)**AB(1, 2) == AB(1, 4)
    assert 2**AB(1, 2) == AB(2, 4)
    assert AB(1, 2)/AB(-1, 1) == AB(-oo, oo)
    assert AB(0, 1)**oo == AB(0, oo)
    assert AB(-2, -1)**-1 == AB(-1, Rational(-1, 2))
    assert str(AB(-1, 1)/0) == 'zoo*AccumBounds(-1, 1)'
    assert str(AB(1, 2)*x) == 'x*AccumBounds(1, 2)'
    assert str(x - AB(1, 2)) == 'x + AccumBounds(-2, -1)'
    assert spy.calls >= 10
    assert abs(AB(-3, 1)) == AB(0, 3)


def test_constructor_errors():
    with pytest.raises(ValueError, match='Lower limit should be smaller'):
        AB(2, 1)
    with pytest.raises(ValueError, match='Only real AccumulationBounds'):
        AB(I, 1)
    with pytest.raises(ValueError, match='Only real AccumulationBounds'):
        AB(x, x + 1)
    assert AB(1, 1) == 1


def test_comparisons_and_sets():
    assert (AB(1, 2) < 3) is symox.true or bool(AB(1, 2) < 3) is True
    assert bool(AB(1, 2) > AB(3, 4)) is False
    held = AB(1, 2) < Rational(3, 2)
    assert str(held) == 'AccumBounds(1, 2) < 3/2'
    assert AB(1, 2).intersection(AB(0, 3)) == AB(1, 2)
    assert AB(1, 3).intersection(FiniteSet(0, 2, 3)) == FiniteSet(2, 3)
    assert AB(1, 2).union(AB(2, 5)) == AB(1, 5)
    assert Rational(3, 2) in AB(1, 2)
    assert 3 not in AB(1, 2)
    assert AB(1, 2).__contains__(Rational(3, 2)) is symox.true
    assert AB(1, 2).__contains__(oo) is False
    with pytest.raises(TypeError, match='input failed to evaluate'):
        Symbol('r', real=True) in AB(1, 2)


def test_function_bounds_are_native(monkeypatch):
    spy = _Spy(monkeypatch, '_accum_evaluate')
    assert sin(AB(0, pi)) == AB(0, 1)
    assert exp(AB(0, 1)) == AB(1, E)
    assert log(AB(1, 2)) == AB(0, log(2))
    assert log(AB(0, 2)) == AB(-oo, log(2))
    # the reference log.eval: a lower end that is not positive or zero
    assert log(AB(-1, 1)) is symox.nan
    assert log(AB(-2, -1)) is symox.nan
    assert cos(oo) == AB(-1, 1)
    assert tan(AB(0, pi)) == AB(-oo, oo)
    assert symox.atan(zoo) == AB(-pi/2, pi/2)
    assert spy.calls >= 6
