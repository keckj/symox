"""Gamma, beta, harmonic and combinatorial rewrites match SymPy 1.14.

Expected strings were produced by SymPy 1.14 for the same calls. They cover
the rule conditions, hints, pattern selection and the bottom-up traversal of
SymPy's ``Basic._rewrite`` through ordinary and structural nodes.
"""
import pytest

import symox as sx

x, y, t = sx.symbols('x y t')
n = sx.Symbol('n', integer=True, nonnegative=True)
k = sx.Symbol('k', integer=True)
m = sx.Symbol('m', integer=True)
p = sx.Symbol('p', positive=True)

RF_GAMMA = 'Piecewise((gamma(x + y)/gamma(x), x > 0), ((-1)**y*gamma(1 - x)/gamma(-x - y + 1), True))'

CASES = [
    (lambda: sx.lowergamma(x, y), 'uppergamma', 'gamma(x) - uppergamma(x, y)'),
    (lambda: sx.lowergamma(x, y), 'expint', '-y**x*expint(1 - x, y) + gamma(x)'),
    (lambda: sx.lowergamma(-2, y), 'expint', 'lowergamma(-2, y)'),
    (lambda: sx.uppergamma(x, y), 'lowergamma', 'gamma(x) - lowergamma(x, y)'),
    (lambda: sx.uppergamma(x, y), 'tractable', 'exp(loggamma(x)) - lowergamma(x, y)'),
    (lambda: sx.uppergamma(x, y), 'expint', 'y**x*expint(1 - x, y)'),
    (lambda: sx.polygamma(2, x), 'zeta', '-2*zeta(3, x)'),
    (lambda: sx.polygamma(3, y), 'harmonic', '-6*harmonic(y - 1, 4) + pi**4/15'),
    (lambda: sx.polygamma(x, y), 'zeta', 'polygamma(x, y)'),
    (lambda: sx.polygamma(x, y), 'harmonic', 'polygamma(x, y)'),
    (lambda: sx.digamma(x), 'harmonic', 'harmonic(x - 1) - EulerGamma'),
    (lambda: sx.trigamma(x), 'zeta', 'zeta(2, x)'),
    (lambda: sx.trigamma(x), 'harmonic', '-harmonic(x - 1, 2) + pi**2/6'),
    (lambda: sx.loggamma(x), 'intractable', 'log(gamma(x))'),
    (lambda: sx.factorial(n), 'gamma', 'gamma(n + 1)'),
    (lambda: sx.factorial(n), 'Product', 'Product(_i, (_i, 1, n))'),
    (lambda: sx.factorial(x), 'Product', 'factorial(x)'),
    (lambda: sx.subfactorial(x), 'uppergamma', 'exp(-1)*uppergamma(x + 1, -1)'),
    (lambda: sx.subfactorial(x), 'gamma',
     '((-1)**(x + 1)*exp(-I*pi*x)*lowergamma(x + 1, -1) + gamma(x + 1))*exp(-1)'),
    (lambda: sx.factorial2(x), 'gamma',
     '2**(x/2)*Piecewise((1, Eq(Mod(x, 2), 0)), (sqrt(2)/sqrt(pi), Eq(Mod(x, 2), 1)))*gamma(x/2 + 1)'),
    (lambda: sx.rf(x, y), 'gamma', RF_GAMMA),
    (lambda: sx.rf(x, k), 'FallingFactorial', 'FallingFactorial(k + x - 1, k)'),
    (lambda: sx.rf(x, k), 'binomial', 'binomial(k + x - 1, k)*factorial(k)'),
    (lambda: sx.rf(x, k), 'factorial', 'RisingFactorial(x, k)'),
    (lambda: sx.rf(m, k), 'factorial',
     'Piecewise((factorial(k + m - 1)/factorial(m - 1), m > 0), ((-1)**k*factorial(-m)/factorial(-k - m), True))'),
    (lambda: sx.ff(n, k), 'gamma', 'gamma(n + 1)/gamma(-k + n + 1)'),
    (lambda: sx.ff(n, k), 'RisingFactorial', 'RisingFactorial(-k + n + 1, k)'),
    (lambda: sx.ff(m, k), 'factorial',
     'Piecewise((factorial(m)/factorial(-k + m), m >= 0), ((-1)**k*factorial(k - m - 1)/factorial(-m - 1), True))'),
    (lambda: sx.binomial(x, k), 'FallingFactorial', 'FallingFactorial(x, k)/factorial(k)'),
    (lambda: sx.binomial(x, y), 'FallingFactorial', 'binomial(x, y)'),
    (lambda: sx.binomial(x, y), 'tractable',
     'exp(loggamma(x + 1))*exp(-loggamma(y + 1))*exp(-loggamma(x - y + 1))'),
    (lambda: sx.beta(x, y), 'gamma', 'gamma(x)*gamma(y)/gamma(x + y)'),
    (lambda: sx.beta(t, y), 'Integral', 'Integral(_t0**(t - 1)*(1 - _t0)**(y - 1), (_t0, 0, 1))'),
    (lambda: sx.harmonic(x), 'digamma', 'polygamma(0, x + 1) + EulerGamma'),
    (lambda: sx.harmonic(x, 2), 'Sum', 'Sum(_k**(-2), (_k, 1, x))'),
    (lambda: sx.harmonic(x, 2), 'zeta', '-zeta(2, x + 1) + pi**2/6'),
    (lambda: sx.harmonic(x, 2), 'tractable', '-polygamma(1, x + 1) + pi**2/6'),
    # Bottom-up traversal through ordinary and structural nodes.
    (lambda: sx.rf(x, y) + sx.gamma(x), 'factorial', 'RisingFactorial(x, y) + factorial(x - 1)'),
    (lambda: sx.rf(x, y) + sx.gamma(x), 'tractable',
     'Piecewise((exp(-loggamma(x))*exp(loggamma(x + y)), x > 0), '
     '((-1)**y*exp(loggamma(1 - x))*exp(-loggamma(-x - y + 1)), True)) + exp(loggamma(x))'),
    (lambda: sx.sin(sx.rf(x, y)), 'gamma', 'sin(%s)' % RF_GAMMA),
    (lambda: sx.Piecewise((sx.ff(x, y), x > 0), (sx.factorial(x), True)), 'gamma',
     'Piecewise((Piecewise((gamma(x + 1)/gamma(x - y + 1), x >= 0), '
     '((-1)**y*gamma(-x + y)/gamma(-x), True)), x > 0), (gamma(x + 1), True))'),
    (lambda: sx.Integral(sx.gamma(t), (t, 0, x)), 'factorial', 'Integral(factorial(t - 1), (t, 0, x))'),
    (lambda: sx.Eq(sx.rf(x, y), sx.gamma(x)), 'gamma', 'Eq(%s, gamma(x))' % RF_GAMMA),
    (lambda: sx.rf(x, y)*sx.sin(x) + sx.beta(x, y), 'gamma',
     '%s*sin(x) + gamma(x)*gamma(y)/gamma(x + y)' % RF_GAMMA),
]


@pytest.mark.parametrize('build, rule, expected', CASES)
def test_rewrite_matches_sympy(build, rule, expected):
    assert str(build().rewrite(rule)) == expected


def test_rewrite_hints_patterns_and_depth():
    assert str(sx.rf(x, y).rewrite('gamma', piecewise=False)) == 'gamma(x + y)/gamma(x)'
    assert str(sx.rf(-p, y).rewrite('gamma', piecewise=False)) == '(-1)**y*gamma(p + 1)/gamma(p - y + 1)'
    assert str(sx.ff(x, y).rewrite(sx.gamma, piecewise=False)) == 'gamma(x + 1)/gamma(x - y + 1)'
    assert str(sx.rf(x, y).rewrite(sx.gamma)) == RF_GAMMA
    # Pattern classes absent from the expression leave it unchanged.
    assert sx.rf(x, y).rewrite(sx.factorial, sx.gamma) == sx.rf(x, y)
    assert str((sx.rf(x, y)+sx.factorial(n)).rewrite(sx.factorial, sx.gamma)) == \
        'RisingFactorial(x, y) + gamma(n + 1)'
    assert str((sx.rf(x, y)+sx.factorial(n)).rewrite([sx.rf, sx.sin], sx.gamma)) == \
        RF_GAMMA+' + factorial(n)'
    assert str(sx.rf(x, sx.rf(x, y)).rewrite('gamma', deep=False)) == (
        'Piecewise((gamma(x + RisingFactorial(x, y))/gamma(x), x > 0), '
        '((-1)**RisingFactorial(x, y)*gamma(1 - x)/gamma(-x - RisingFactorial(x, y) + 1), True))')
    assert str(sx.rf(x, k).rewrite('tractable', limitvar=k)) == 'exp(loggamma(k + x))/gamma(x)'
    assert str(sx.ff(x, -k).rewrite('tractable', limitvar=k)) == 'exp(-loggamma(k + x + 1))*gamma(x + 1)'


@pytest.mark.parametrize('call,expected', [
    (lambda: sx.fibonacci(x).rewrite(sx.sqrt), 'sqrt(5)*(-(1 - sqrt(5))**x + (1 + sqrt(5))**x)/(5*2**x)'),
    (lambda: sx.fibonacci(x).rewrite('sqrt'), 'sqrt(5)*(-(1 - sqrt(5))**x + (1 + sqrt(5))**x)/(5*2**x)'),
    (lambda: sx.fibonacci(x).rewrite(sx.GoldenRatio), '(GoldenRatio**x - 1/(-GoldenRatio)**x)/(-1 + 2*GoldenRatio)'),
    (lambda: sx.fibonacci(x).rewrite('tractable'), 'sqrt(5)*(GoldenRatio**x - cos(pi*x)/GoldenRatio**x)/5'),
    (lambda: sx.lucas(x).rewrite(sx.sqrt), '((1 - sqrt(5))**x + (1 + sqrt(5))**x)/2**x'),
    (lambda: sx.tribonacci(x).rewrite(sx.TribonacciConstant),
     'floor(3*TribonacciConstant**x*(102*sqrt(33) + 586)**(1/3)/(-2*(102*sqrt(33) + 586)**(1/3) + 4 + (102*sqrt(33) + 586)**(2/3)) + 1/2)'),
    (lambda: (sx.fibonacci(x) + sx.lucas(x)).rewrite(sx.sqrt),
     'sqrt(5)*(-(1 - sqrt(5))**x + (1 + sqrt(5))**x)/(5*2**x) + ((1 - sqrt(5))**x + (1 + sqrt(5))**x)/2**x'),
    (lambda: sx.bernoulli(x).rewrite(sx.zeta), 'Piecewise((1, Eq(x, 0)), (-x*zeta(1 - x), True))'),
    (lambda: sx.bernoulli(x, y).rewrite(sx.zeta), 'Piecewise((1, Eq(x, 0)), (-x*zeta(1 - x, y), True))'),
    (lambda: sx.euler(x, y).rewrite(sx.genocchi),
     'Piecewise((-polygamma(0, y/2) + polygamma(0, y/2 + 1/2), Eq(x, -1)), (-genocchi(x + 1, y)/(x + 1), True))'),
    (lambda: sx.andre(x).rewrite(sx.zeta),
     '2*(2*pi)**(-x - 1)*(-cos(pi*x)*zeta(x + 1, 3/4) + zeta(x + 1, 1/4))*gamma(x + 1)'),
    (lambda: sx.andre(x).rewrite(sx.polylog), 'I**(x + 1)*polylog(-x, -I) + (-I)**(x + 1)*polylog(-x, I)'),
    (lambda: sx.bernoulli(x, 1), 'bernoulli(x)'),
    (lambda: sx.mathematica_code(sx.fibonacci(x)), '(1/5)*2^(-x)*5^(1/2)*(-(1 - 5^(1/2))^x + (1 + 5^(1/2))^x)'),
])
def test_fibonacci_family_and_number_rewrites(call, expected):
    assert str(call()) == expected


def test_tribonacci_sqrt_rewrite_is_the_closed_form():
    form = sx.tribonacci(x).rewrite(sx.sqrt)
    assert form.has(sx.I) and not form.has(sx.tribonacci)
    for value in (1, 2, 7, 10):
        assert abs(form.subs(x, value).evalf(30) - sx.tribonacci(value)) < sx.Float('1e-25')


def test_genocchi_rewrites():
    assert str(sx.genocchi(x).rewrite(sx.dirichlet_eta)) == '-2*x*dirichlet_eta(1 - x)'
    # SymPy 1.14: -2*2**x*bernoulli(x, 1.0) + 2*bernoulli(x); the default
    # x = 1 is kept exact here
    assert str(sx.genocchi(x).rewrite(sx.bernoulli)) == '-2*2**x*bernoulli(x) + 2*bernoulli(x)'
    assert str(sx.genocchi(n).rewrite(sx.bernoulli)) == '(2 - 2*2**n)*bernoulli(n)'
