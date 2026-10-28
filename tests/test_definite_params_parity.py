"""Definite integrals with parameters follow the reference only where its
value holds for every parameter value the assumptions allow."""
import pytest
import symox as sx
import sympy as sp

ASSUMPTIONS = {
    'none': {},
    'pos': {'positive': True},
    'real': {'real': True},
    'cplx': {'complex': True},
    'nz': {'nonzero': True},
}


def evaluate(module, assumption, call):
    names = {k: getattr(module, k) for k in ('integrate', 'limit', 'log', 'oo', 'Heaviside', 'DiracDelta', 'Integral')}
    names['x'] = module.Symbol('x')
    names['a'] = module.Symbol('a', **ASSUMPTIONS[assumption])
    names['b'] = module.Symbol('b', **ASSUMPTIONS[assumption])
    return eval(call, names)


# cases where the reference value is right for all allowed parameters
AGREED = [
    ('pos', 'integrate(1/(x+a), (x, 0, 1))'),
    ('pos', 'integrate(1/(1+a*x), (x, 0, 1))'),
    ('pos', 'integrate(1/x, (x, a, 1))'),
    ('pos', 'integrate(x**a, (x, 1, oo))'),
    ('none', 'integrate(log(a*x), (x, 0, 1))'),
    ('real', 'integrate(log(a*x), (x, 0, 1))'),
    ('cplx', 'integrate(log(a*x), (x, 0, 1))'),
    ('pos', 'integrate(log(a*x), (x, 0, 1))'),
    ('nz', 'integrate(log(a*x), (x, 0, 1))'),
    ('none', 'integrate(log(a*x), x)'),
    ('nz', 'integrate(Heaviside(x-a), (x, 0, 1))'),
    ('pos', 'integrate(Heaviside(x-a), (x, 0, 1))'),
    ('none', 'integrate(DiracDelta(x-a), (x, 0, 1))'),
    ('pos', 'integrate(DiracDelta(x-a), (x, 0, 1))'),
    ('pos', 'limit(x**a, x, oo)'),
    ('pos', 'limit(x**(-a), x, oo)'),
    ('pos', 'limit(x**(a+1)/(a+1), x, oo)'),
]


@pytest.mark.parametrize('assumption,call', AGREED)
def test_agreed_reference_values(assumption, call):
    actual = evaluate(sx, assumption, call).to_sympy()
    expected = evaluate(sp, assumption, call)
    assert not actual.has(sp.Integral) or expected.has(sp.Integral)
    assert actual == expected or sp.simplify(actual - expected) == 0


# The reference returns these values without a condition, although the
# integrand has a pole on the interval (or a divergence) for some allowed
# parameter values, e.g. -log(a) for 1/x on (a, 1) with a = -1. They stay
# unevaluated (see docs/algorithmic-parity.md).
HELD = [
    (assumption, call)
    for assumption in ('none', 'real', 'cplx', 'nz')
    for call in (
        'integrate(1/(x+a), (x, 0, 1))',
        'integrate(1/(x-a), (x, 0, 1))',
        'integrate(1/(1+a*x), (x, 0, 1))',
        'integrate(1/x, (x, a, 1))',
    )
] + [('pos', 'integrate(1/(x-a), (x, 0, 1))')]


@pytest.mark.parametrize('assumption,call', HELD)
def test_values_wrong_for_some_parameters_stay_held(assumption, call):
    actual = evaluate(sx, assumption, call)
    if not (assumption=='nz' and call=='integrate(1/(1+a*x), (x, 0, 1))') and call in ('integrate(1/(x+a), (x, 0, 1))',
                'integrate(1/(x-a), (x, 0, 1))',
                'integrate(1/(1+a*x), (x, 0, 1))',
                'integrate(1/x, (x, a, 1))'):
        # A proven pole-free branch is now available. The full original
        # Integral remains the default on its excluded parameter domain.
        assert actual.is_Piecewise
        fallback,condition=actual.args[-1]
        assert condition==sx.true
        assert fallback==evaluate(sx,assumption,call.replace('integrate(', 'Integral('))
    else:
        assert isinstance(actual.to_sympy(), sp.Integral)


@pytest.mark.parametrize('integrand,limits', [
    ('1/(x+a)', '(x, 2, 3)'),
    ('1/(x-a)', '(x, -2, -1)'),
    ('1/(x+a)**2', '(x, 0, 1)'),
    ('1/(x**2+a)', '(x, 0, 1)'),
])
def test_poles_outside_the_interval_by_assumption(integrand, limits):
    actual = sx.integrate(*evaluate(sx, 'pos', f'({integrand}, {limits})')).to_sympy()
    expected = sp.integrate(*evaluate(sp, 'pos', f'({integrand}, {limits})'))
    assert not actual.has(sp.Integral)
    a = sp.Symbol('a', positive=True)
    for value in (sp.Rational(7, 10), sp.Rational(23, 10)):
        assert abs(complex((actual - expected).subs(a, value).evalf(30))) < 1e-20


def test_pole_inside_by_assumption_diverges():
    a = sx.Symbol('a', positive=True)
    x = sx.Symbol('x')
    # the pole -a lies in (-a-1, 0): one-sided infinite contributions
    assert sx.integrate(1/(x+a)**2, (x, -a-1, 0)).to_sympy() == sp.oo
