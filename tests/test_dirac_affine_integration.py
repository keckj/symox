"""Exact compact-support distribution rules, their domains and held fallbacks."""
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize('slope', [1, -1, 2, -2])
@pytest.mark.parametrize('root', [0, 1, -2])
@pytest.mark.parametrize('weight', ['one', 'linear', 'quadratic', 'exponential', 'cosine', 'rational'])
def test_affine_delta_with_analytic_multiplier_matches_reference(slope, root, weight):
    def run(module):
        t = module.Symbol('t', real=True)
        rate = module.Symbol('s')
        multiplier = dict(one=module.Integer(1), linear=t, quadratic=t*t,
            exponential=module.exp(-rate*t), cosine=module.cos(t), rational=1/(1+t*t))[weight]
        return module.integrate(multiplier*module.DiracDelta(slope*(t-root)), (t, -5, 5))
    actual, expected = run(sx), run(sp)
    assert sp.simplify(actual.to_sympy()-expected) == 0
    assert not actual.has(sx.Integral)


@pytest.mark.parametrize('bounds', [(-2, 2), (2, -2), (1, 2)])
def test_support_and_endpoint_half_weight_are_preserved(bounds):
    t = sx.Symbol('t', real=True)
    x = sp.Symbol('t', real=True)
    actual = sx.integrate(sx.exp(t)*sx.DiracDelta(t), (t, *bounds))
    expected = sp.integrate(sp.exp(x)*sp.DiracDelta(x), (x, *bounds))
    assert actual.to_sympy() == expected


@pytest.mark.parametrize('bounds', [(0, 2), (-2, 0)])
def test_finite_origin_endpoints_use_symmetric_regularization_and_are_not_reference_matches(bounds):
    t = sx.Symbol('t', real=True)
    x = sp.Symbol('t', real=True)
    actual = sx.integrate(sx.exp(t)*sx.DiracDelta(t), (t, *bounds))
    assert actual == sx.Rational(1, 2)
    # SymPy 1.14 returns zero here, but half weight after translation or at
    # an infinite endpoint. Record the difference rather than hiding it.
    assert sp.integrate(sp.exp(x)*sp.DiracDelta(x), (x, *bounds)) == 0
    shifted = sx.integrate(sx.exp(t-1)*sx.DiracDelta(t-1), (t, bounds[0]+1, bounds[1]+1))
    assert shifted == actual
    epsilon = sp.Symbol('epsilon', positive=True)
    # Independent normalized Gaussian regularization: integral_0^2 delta_eps.
    assert sp.limit(sp.erf(2/epsilon)/2, epsilon, 0, dir='+') == sp.Rational(1, 2)


def test_unilateral_limits_need_no_exponential_convergence_assumption():
    t = sx.Symbol('t', real=True)
    rate = sx.Symbol('s')
    lower = sx.Symbol('t0', negative=True)
    expression = sx.exp(-rate*t)*sx.DiracDelta(t)
    assert sx.integrate(expression, (t, lower, sx.oo)) == 1
    assert sx.integrate(expression, (t, 0, sx.oo)) == sx.Rational(1, 2)


@pytest.mark.parametrize('value', ['zero', 'one', 'minus_one', 'infinity', 'minus_infinity',
    'positive', 'negative', 'extended_positive', 'extended_negative'])
@pytest.mark.parametrize('at_zero', [0, 1, 'parameter', 'nan'])
def test_heaviside_custom_zero_value_and_extended_real_signs(value, at_zero):
    def run(module):
        argument = dict(zero=module.Integer(0), one=module.Integer(1), minus_one=module.Integer(-1),
            infinity=module.oo, minus_infinity=-module.oo,
            positive=module.Symbol('p', positive=True), negative=module.Symbol('n', negative=True),
            extended_positive=module.Symbol('ep', extended_positive=True),
            extended_negative=module.Symbol('en', extended_negative=True))[value]
        choice = module.Symbol('h0') if at_zero == 'parameter' else module.nan if at_zero == 'nan' else at_zero
        return module.Heaviside(argument, choice)
    assert run(sx).to_sympy() == run(sp)


@pytest.mark.parametrize('name', ['Heaviside', 'DiracDelta'])
def test_proved_complex_arguments_raise_and_explicitly_held_arguments_remain_held(name):
    for module in (sx, sp):
        function = getattr(module, name)
        with pytest.raises(ValueError):
            function(1+module.I)
        assert function(1+module.I, evaluate=False).args[0] == 1+module.I


@pytest.mark.parametrize('order', [-1, .5, 1.0, 'symbol'])
def test_dirac_derivative_order_requires_nonnegative_integer_literal(order):
    for module in (sx, sp):
        t = module.Symbol('t', real=True)
        n = module.Symbol('n', integer=True, nonnegative=True) if order == 'symbol' else order
        with pytest.raises(ValueError):
            module.DiracDelta(t, n)
        assert len(module.DiracDelta(t, n, evaluate=False).args) == 2


def test_singular_and_multiple_distribution_products_are_not_assigned_values():
    t = sx.Symbol('t', real=True)
    for expression in [sx.DiracDelta(t)/t,
                       sx.DiracDelta(t)**2, sx.DiracDelta(t*t-1)]:
        assert sx.integrate(expression, t).has(sx.Integral)


def test_existing_substitution_rule_for_step_delta_products_is_preserved():
    t = sx.Symbol('t', real=True)
    # This uses the established u=Heaviside(t) substitution, not the new
    # analytic multiplier rule. It must remain available for inexact input.
    exact = sx.integrate(sx.Heaviside(t)*sx.DiracDelta(t), (t, -1, 1))
    inexact = sx.integrate(sx.Float(2)*sx.Heaviside(t)*sx.DiracDelta(t), (t, -1, 1))
    assert exact == sx.Rational(1, 2)
    assert inexact == sx.Float(1)


def test_inexact_dependent_delta_retains_held_fallback_without_underflow():
    t = sx.Symbol('t', real=True)
    expression = sx.exp(t)*sx.DiracDelta(sx.Float('1e300')*t-sx.Float('1e-300'))
    assert sx.integrate(expression, t).has(sx.Integral)
    # A delta independent of the variable keeps its established constant path.
    result = sx.integrate(sx.Float(2)*sx.DiracDelta(0), t)
    assert not result.has(sx.Integral)


def test_large_numeric_powers_are_preflighted_but_cheap_large_powers_remain_supported():
    script = '''
import symox as s
t=s.Symbol('t',real=True)
power=s.Integer(2)**80
held=s.integrate(t**power*s.DiracDelta(t-2),t)
assert held.has(s.Integral)
cheap=s.integrate(t**power*s.DiracDelta(t-1),t)
assert cheap==s.Heaviside(t-1)
'''
    subprocess.run([sys.executable, '-c', script], check=True, timeout=5)
