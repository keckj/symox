"""Integration results that print exactly like the reference.

The expected strings were produced by the reference SymPy 1.14 on the same
inputs: Piecewise integrands (indefinite, definite and improper), the
degenerate cases of polynomial times exponential or trigonometric
antiderivatives, Meijer G indefinite integration, the branch factors of the
incomplete gamma and exponential integrals and logarithm expansion.
"""
import pytest
import symox
from symox import *
from symox.integrals.meijerint import meijerint_indefinite

x, y = symbols('x y')
a, b = symbols('a b')
r = Symbol('r', real=True)
p, q = symbols('p q', positive=True)

CASES = [   ('integrate(Piecewise((x, x < 1), (x**2, True)), (x, 0, 2))', '17/6'),
    (   'integrate(Piecewise((x, x < 1), (x**2, True)), x)',
        'Piecewise((x**2/2, x < 1), (x**3/3 + 1/6, True))'),
    (   'integrate(Piecewise((x**2, x < y), (5*x**4, True)), (x, 0, 1))',
        '-Min(1, Max(0, y))**5 + Min(1, Max(0, y))**3/3 + 1'),
    (   'integrate(Piecewise((x**3/3, Abs(x) < 1), (x**5, True)).diff(x), x)',
        'Piecewise((x**5, x <= -1), (x**3/3 - 2/3, x < 1), (x**5 - 4/3, True))'),
    (   'integrate(Piecewise((0, x < 0), (1, x < 1), (2, True)), x)',
        'Piecewise((0, x < 0), (x, x < 1), (2*x - 1, True))'),
    (   'Piecewise((0, x < 0), (1, x < 1), (2, True)).piecewise_integrate(x)',
        'Piecewise((0, x < 0), (x, x < 1), (2*x, True))'),
    ('integrate(Abs(x), (x, -1, 2))', '5/2'),
    ('integrate(Abs(x), x)', 'Integral(Abs(x), x)'),
    ('integrate(Abs(r), r)', 'Piecewise((-r**2/2, r <= 0), (r**2/2, True))'),
    ('integrate(Max(x, 1), (x, 0, 2))', '5/2'),
    ('integrate(Min(x, 1)*x, (x, 0, 3))', '13/3'),
    ('integrate(sign(x)*x, (x, -1, 1))', '1'),
    ('integrate(Abs(x - 1), (x, 0, 3))', '5/2'),
    ('integrate(Piecewise((1, x > 0), (0, True)), (x, -oo, 1))', '1'),
    ('integrate(Piecewise((exp(-x), x > 0), (0, True)), (x, -oo, oo))', '1'),
    ('integrate(Piecewise((1, x < p), (0, True)), (x, 0, 2*p))', 'p'),
    ('integrate(Piecewise((1, (x > 0) & (x < 1)), (0, True)), (x, -2, 2))', '1'),
    ('integrate(Piecewise((1, (x < 0) | (x > 1)), (0, True)), (x, -2, 2))', '3'),
    ('integrate(Piecewise((sin(x), x < pi), (0, True)), (x, 0, 2*pi))', '2'),
    ('integrate(Piecewise((1, Ne(x, 0)), (x, True)), (x, -1, 1))', '2'),
    ('integrate(Piecewise((x, x < 1), (1/x, True)), (x, 0, oo))', 'oo'),
    ('integrate(Piecewise((x, x < 1), (1/x**2, True)), (x, 0, oo))', '3/2'),
    (   'integrate(Piecewise((x, y > 0), (x**2, True)), (x, 0, 1))',
        'Piecewise((1/2, (y > 0) & (y < oo)), (1/3, True))'),
    ('integrate(Piecewise((1, x**2 < 1), (0, True)), (x, -2, 2))', '2'),
    ('integrate(Piecewise((x, x < 1), (x**2, True)), (x, 2, 0))', '-17/6'),
    (   'integrate(Piecewise((x, x < 1), (x**2, True)), (x, 0, y))',
        'Piecewise((y**2/2, y < 0), (y**3/3 - Min(1, y)**3/3 + Min(1, y)**2/2, True))'),
    ('integrate(Piecewise((x, x < 1)), x)', 'Piecewise((x**2/2, x < 1), (nan, True))'),
    ('integrate(Max(0, 1 - Abs(x)), (x, -2, 2))', '1'),
    ('integrate(Piecewise((x, x < 1), (2, True))**2, (x, 0, 3))', '25/3'),
    ('Integral(Piecewise((x, x < 1), (2, True)), (x, 0, 3)).doit()', '9/2'),
    ('integrate(Piecewise((cos(x), x < 0), (exp(-x), True)), (x, -pi, oo))', '1'),
    (   'integrate(x*exp(a*x), (x, 0, y))',
        'Piecewise(((a*y - 1)*exp(a*y)/a**2 + a**(-2), (a > -oo) & (a < oo) & Ne(a, 0)), (y**2/2, '
        'True))'),
    (   'integrate(x*exp(a*x), x)',
        'Piecewise(((a*x - 1)*exp(a*x)/a**2, Ne(a**2, 0)), (x**2/2, True))'),
    (   'integrate(x**2*exp(a*x), x)',
        'Piecewise(((a**2*x**2 - 2*a*x + 2)*exp(a*x)/a**3, Ne(a**3, 0)), (x**3/3, True))'),
    (   'integrate(x*exp(2*a*x), x)',
        'Piecewise(((2*a*x - 1)*exp(2*a*x)/(4*a**2), Ne(a**2, 0)), (x**2/2, True))'),
    (   'integrate(x*exp(I*a*x), x)',
        'Piecewise(((-I*a*x + 1)*exp(I*a*x)/a**2, Ne(a**2, 0)), (x**2/2, True))'),
    (   'integrate(x*exp(a*b*x), x)',
        'Piecewise(((a*b*x - 1)*exp(a*b*x)/(a**2*b**2), Ne(a**2*b**2, 0)), (x**2/2, True))'),
    (   'integrate(x*exp((a + b)*x), x)',
        'Piecewise(((a*x + b*x - 1)*exp(x*(a + b))/(a**2 + 2*a*b + b**2), Ne(a**2 + 2*a*b + b**2, '
        '0)), (x**2/2, True))'),
    (   'integrate(x**2*exp(-a*x/3), x)',
        'Piecewise(((-3*a**2*x**2 - 18*a*x - 54)*exp(-a*x/3)/a**3, Ne(a**3, 0)), (x**3/3, True))'),
    ('integrate(x*exp(2*x), x)', '(2*x - 1)*exp(2*x)/4'),
    ('integrate(x**3*exp(p*x), x)', '(p**3*x**3 - 3*p**2*x**2 + 6*p*x - 6)*exp(p*x)/p**4'),
    ('integrate(x*exp(sqrt(2)*x), x)', '(sqrt(2)*x - 1)*exp(sqrt(2)*x)/2'),
    ('integrate(x*sin(a*x), x)', 'Piecewise((-x*cos(a*x)/a + sin(a*x)/a**2, Ne(a, 0)), (0, True))'),
    (   'integrate(x*cos(a*x), x)',
        'Piecewise((x*sin(a*x)/a + cos(a*x)/a**2, Ne(a, 0)), (x**2/2, True))'),
    (   'integrate(x*sin(2*a*x), x)',
        'Piecewise((-x*cos(2*a*x)/(2*a) + sin(2*a*x)/(4*a**2), Ne(a, 0)), (0, True))'),
    (   'integrate(x*cos(a*b*x), x)',
        'Piecewise((x**2/2, Eq(a, 0) | Eq(b, 0)), (x*sin(a*b*x)/(a*b) + cos(a*b*x)/(a**2*b**2), '
        'True))'),
    (   'integrate(x**2*cos(a*x), x)',
        'Piecewise((x**2*sin(a*x)/a + 2*x*cos(a*x)/a**2 - 2*sin(a*x)/a**3, Ne(a, 0)), (x**3/3, '
        'True))'),
    ('integrate(exp(a*x), x)', 'Piecewise((exp(a*x)/a, Ne(a, 0)), (x, True))'),
    ('integrate(sin(a*x), x)', 'Piecewise((-cos(a*x)/a, Ne(a, 0)), (0, True))'),
    (   'integrate(exp(a*x)*sin(a*x), x)',
        'Piecewise((exp(a*x)*sin(a*x)/(2*a) - exp(a*x)*cos(a*x)/(2*a), Ne(a, 0)), (0, True))'),
    ('meijerint_indefinite(x*besselj(0, x), x)', 'x*besselj(1, x)'),
    ('meijerint_indefinite(exp(-x**2), x)', 'sqrt(pi)*erf(x)/2'),
    ('meijerint_indefinite(1/(x - 1)**2, x)', '-x/(x - 1)'),
    ('meijerint_indefinite(besselj(1, x)**2/x, x)', '-besselj(0, x)**2/2 - besselj(1, x)**2/2'),
    ('meijerint_indefinite(sin(x)/x, x)', 'Si(x)'),
    ('meijerint_indefinite(x**2*exp(x), x)', 'x**2*exp(x) - 2*x*exp(x) + 2*exp(x)'),
    ('meijerint_indefinite(x*exp(-x**2), x)', '-exp(-x**2)/2'),
    ('meijerint_indefinite(1/(1 + x**2), x)', 'atan(x)'),
    ('meijerint_indefinite(sin(x)**2, x)', 'x/2 - sin(2*x)/4'),
    ('meijerint_indefinite(erf(x), x)', 'x*erf(x) + exp(-x**2)/sqrt(pi)'),
    ('integrate(sin(x)/x, x, meijerg=True)', 'Si(x)'),
    ('integrate(x**2, (x, 0, 1), meijerg=True)', '1/3'),
    ('Integral(exp(-x**2), (x, 0, oo)).doit(meijerg=True)', 'sqrt(pi)/2'),
    ('integrate(exp(-x**2)*x**3, x)', '(-x**2 - 1)*exp(-x**2)/2'),
    ('lowergamma(a, x*exp_polar(2*I*pi))', 'exp(2*I*pi*a)*lowergamma(a, x)'),
    (   'uppergamma(a, x*exp_polar(2*I*pi))',
        '(1 - exp(2*I*pi*a))*gamma(a) + exp(2*I*pi*a)*uppergamma(a, x)'),
    ('uppergamma(3, x*exp_polar(2*I*pi))', '2*(x**2/2 + x + 1)*exp(-x)'),
    ('Ei(x*exp_polar(2*I*pi))', 'Ei(x) + 2*I*pi'),
    (   'expint(a, x*exp_polar(2*I*pi))',
        'x**(a - 1)*(exp(2*I*pi*a) - 1)*gamma(1 - a) + expint(a, x)'),
    ('expint(2, x*exp_polar(2*I*pi))', '2*I*pi*x + expint(2, x)'),
    ('expand_log(log(p**2*q))', '2*log(p) + log(q)'),
    ('expand_log(log(p*q)**2)', '(log(p) + log(q))**2'),
    ('expand_log(log(p/q) + log(p**3), force=False)', '4*log(p) - log(q)'),
    ('expand_log(log(x*y), force=True)', 'log(x) + log(y)'),
    ('piecewise_fold(x*Piecewise((1, x < 0), (x, True)))', 'Piecewise((x, x < 0), (x**2, True))')]


@pytest.mark.parametrize('source, expected', CASES, ids=[c for c, _ in CASES])
def test_prints_like_the_reference(source, expected):
    assert str(eval(source)) == expected


# exp(a*x)*trig(b*x): the parameter values that zero the denominator of the
# generic antiderivative are kept as separate cases, like the reference
# heurisch_wrapper. The special antiderivatives are checked by derivative
# because their form differs from the reference heurisch.
@pytest.mark.parametrize('integrand, conditions', [
    (lambda: exp(a*x)*sin(b*x), ['Eq(a, 0) & Eq(b, 0)', 'Eq(a, -I*b)', 'Eq(a, I*b)', 'True']),
    (lambda: exp(a*x)*cos(x), ['Eq(a, -I)', 'Eq(a, I)', 'True']),
    (lambda: x*exp(a*x)*sin(x), ['Eq(a, -I)', 'Eq(a, I)', 'True']),
    (lambda: x*exp(a*x)*sin(b*x), ['Eq(a, 0) & Eq(b, 0)', 'Eq(a, -I*b)', 'Eq(a, I*b)', 'True']),
    (lambda: exp(2*a*x)*cos(b*x), ['Eq(a, 0) & Eq(b, 0)', 'Eq(a, -I*b/2)', 'Eq(a, I*b/2)', 'True']),
    (lambda: exp(a*x)*sin(2*x), ['Eq(a, -2*I)', 'Eq(a, 2*I)', 'True']),
])
def test_exp_trig_degenerate_cases(integrand, conditions):
    f = integrand()
    result = integrate(f, x)
    assert isinstance(result, Piecewise)
    assert [str(c) for _, c in result.args] == conditions
    for value, cond in result.args:
        g = f
        if cond != True:
            for eq in (cond.args if isinstance(cond, And) else [cond]):
                g = g.subs(eq.lhs, eq.rhs)
                value = value.subs(eq.lhs, eq.rhs)
        assert_antiderivative(value, g)


def assert_antiderivative(value, f):
    point = {x: Float('0.3', 30), a: Float('0.7', 30), b: Float('1.1', 30)}
    residue = (diff(value, x) - f).subs(point).evalf(30)
    assert abs(complex(residue)) < 1e-20


def test_cyclic_rule_never_divides_by_zero():
    for f in [exp(I*x)*sin(x), exp(-I*x)*cos(x), x*exp(I*x)*sin(x)]:
        result = integrate(f, x)
        assert not result.has(zoo, nan)
        assert_antiderivative(result, f)


def test_quadratic_roots_take_squares_out_of_the_radical():
    assert solve(a**2 + b**2, a) == [-I*b, I*b]
