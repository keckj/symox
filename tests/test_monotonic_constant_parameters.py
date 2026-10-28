"""A parameter's conditional zero set is not a spatial turning point."""
import pytest
import symox as sx


@pytest.mark.parametrize('slope', ['a', 'a**2', '-a**2', 'sin(a)', '0', '3'])
def test_every_finite_real_constant_slope_is_weakly_monotonic(slope):
    x = sx.Symbol('x')
    a, b = sx.symbols('a b', real=True)
    coefficient = sx.sympify(slope, locals={'a': a})
    for domain in (sx.Reals, sx.Interval(-2, 5), sx.Interval.open(0, sx.oo)):
        assert sx.is_monotonic(coefficient*x + b, domain, x) is True


def test_constant_equation_preserves_the_complete_parameter_condition():
    x, a = sx.symbols('x a')
    actual = sx.solveset(a, x)
    assert actual == sx.ConditionSet(x, sx.Eq(a, 0), sx.Complexes)
    assert actual.subs(a, 0) == sx.Complexes
    assert actual.subs(a, sx.I) == sx.EmptySet


def test_held_undefined_constant_is_not_an_affine_monotonicity_proof():
    from symox import _core
    x = sx.Symbol('x')
    coefficient = sx.Heaviside(0, sx.nan, evaluate=False)
    assert sx.Heaviside(*coefficient.args) is sx.nan
    assert not _core._st_constant_real_derivative(coefficient, x)
    assert sx.is_monotonic(coefficient*x, sx.Reals, x) is not True


def test_class_real_fact_does_not_certify_a_held_unknown_value():
    from symox import _core
    x = sx.Symbol('x')
    coefficient = sx.Heaviside(0, sx.Function('f')(1), evaluate=False)
    assert not _core._st_constant_real_derivative(coefficient, x)
    assert sx.is_monotonic(coefficient*x, sx.Reals, x) is not True
