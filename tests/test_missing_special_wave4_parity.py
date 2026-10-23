"""Wave 4 parity: special-function evaluation and rewrites.

Expected values are hard-coded from SymPy 1.14.
"""
import pytest

import symox
from symox.abc import a, b, m, n, s, x, y, z, nu, theta, phi  # noqa: F401

# A star import would also bind symox.test, which pytest would collect.
NS = {k: getattr(symox, k) for k in dir(symox) if not k.startswith('_') and k not in ('test', 'doctest')}
NS.update(a=a, b=b, m=m, n=n, s=s, x=x, y=y, z=z, nu=nu, theta=theta, phi=phi)
globals().update({k: v for k, v in NS.items() if not k.startswith('test')})
from symox.core.function import ArgumentIndexError

CASES = [
    ('erfi(-I*oo)', '-I'),
    ('erfi(I*oo)', 'I'),
    ('erf(-I*oo)', '-oo*I'),
    ('polygamma(0, -oo)', 'oo'),
    ('polygamma(0, I*oo)', 'oo'),
    ('loggamma(I*oo)', 'zoo'),
    ('uppergamma(x, oo)', '0'),
    ('bernoulli(3)', '0'),
    ('bernoulli(1000001)', '0'),
    ("bernoulli(Symbol('k', odd=True, positive=True) + 2)", '0'),
    ('subfactorial(oo)', 'oo'),
    ('lucas(oo)', 'oo'),
    ('fibonacci(oo)', 'oo'),
    ('bell(oo)', 'oo'),
    ('airybi(-oo)', '0'),
    ('airybiprime(-oo)', '0'),
    ('bessely(0, I*oo)', 'oo*I'),
    ('besselk(0, oo)', '0'),
    ('besselk(1, I*oo)', '0'),
    ('conjugate(Ynm(n, m, theta, phi))', '(-1)**(2*m)*exp(-2*I*m*phi)*Ynm(n, m, theta, phi)'),
    ('Ynm_c(n, m, theta, phi)', '(-1)**(2*m)*exp(-2*I*m*phi)*Ynm(n, m, theta, phi)'),
    ('Ynm_c(1, 0, theta, phi)', 'Ynm(1, 0, theta, phi)'),
    ('conjugate(Ynm(2, 1, theta, phi))', 'exp(-2*I*phi)*Ynm(2, 1, theta, phi)'),
    ('SingularityFunction(x, a, nan)', 'nan'),
    ('SingularityFunction(x, nan, 2)', 'nan'),
    ('SingularityFunction(x, a, S(1)/2)', 'SingularityFunction(x, a, 1/2)'),
    ('SingularityFunction(x, a, 3).rewrite(Heaviside)', '(-a + x)**3*Heaviside(-a + x, 1)'),
    ('(1 + SingularityFunction(x, a, -1)).rewrite(DiracDelta)', 'DiracDelta(-a + x) + 1'),
    ('(x*SingularityFunction(x, a, 2)).rewrite(Piecewise)', 'x*Piecewise(((-a + x)**2, a - x <= 0), (0, True))'),
    ('SingularityFunction(x, a, -3).rewrite(Piecewise)', 'Piecewise((oo, Eq(a - x, 0)), (0, True))'),
    ('SingularityFunction(x, a, n).rewrite(Piecewise)', 'SingularityFunction(x, a, n)'),
    ('Heaviside(x).fdiff()', 'DiracDelta(x)'),
    ('DiracDelta(x).fdiff()', 'DiracDelta(x, 1)'),
    ('DiracDelta(x, 2).fdiff()', 'DiracDelta(x, 3)'),
    ('Eq(2, 5, evaluate=False).doit()', 'False'),
    ('Ne(2, 5, evaluate=False).doit()', 'True'),
    ('Lt(2, 5, evaluate=False).doit()', 'True'),
    ('Eq(x, x, evaluate=False).doit()', 'True'),
    ('Si(z).rewrite(expint)', '-I*(-expint(1, z*exp_polar(-I*pi/2))/2 + expint(1, z*exp_polar(I*pi/2))/2) + pi/2'),
    ('Ci(z).rewrite(expint)', '-expint(1, z*exp_polar(-I*pi/2))/2 - expint(1, z*exp_polar(I*pi/2))/2'),
    ('Shi(z).rewrite(expint)', 'expint(1, z)/2 - expint(1, z*exp_polar(I*pi))/2 - I*pi/2'),
    ('Chi(z).rewrite(expint)', '-expint(1, z)/2 - expint(1, z*exp_polar(I*pi))/2 - I*pi/2'),
    ('Si(z).rewrite(sinc)', 'Integral(sinc(_t), (_t, 0, z))'),
    ('Ci(z).rewrite(Ei)', 'Ei(z*exp_polar(-I*pi/2))/2 + Ei(z*exp_polar(I*pi/2))/2'),
    ('Ei(z).rewrite(expint)', '-expint(1, z*exp_polar(I*pi)) - I*pi'),
    ('Ei(z).rewrite(li)', 'li(exp(z))'),
    ('Ei(z).rewrite(Si)', 'Chi(z) + Shi(z)'),
    ('expint(nu, z).rewrite(uppergamma)', 'z**(nu - 1)*uppergamma(1 - nu, z)'),
    ('expint(1, z).rewrite(Ei)', '-Ei(z*exp_polar(I*pi)) + I*pi'),
    ('expint(2, z).rewrite(Ei)', '-z*(-Ei(z*exp_polar(I*pi)) + I*pi) + exp(-z)'),
    ('expint(S(1)/2, z).rewrite(erf)', 'sqrt(pi)*(1 - erf(sqrt(z)))/sqrt(z)'),
    ('li(z).rewrite(Li)', 'Li(z) + li(2)'),
    ('li(z).rewrite(Ei)', 'Ei(log(z))'),
    ('li(z).rewrite(Si)', '-log(I*log(z)) - log(1/log(z))/2 + log(log(z))/2 + Ci(I*log(z)) + Shi(log(z))'),
    ('li(z).rewrite(hyper)', 'log(z)*hyper((1, 1), (2, 2), log(z)) - log(1/log(z))/2 + log(log(z))/2 + EulerGamma'),
    ('Li(z).rewrite(li)', 'li(z) - li(2)'),
    ('erfc(sqrt(z)).rewrite(erf)', '1 - erf(sqrt(z))'),
    ('airyai(z).rewrite(besselj)', 'airyai(z)'),
    ('airybi(z).rewrite(besseli)', 'sqrt(3)*(z*besseli(1/3, 2*z**(3/2)/3)/(z**(3/2))**(1/3) + (z**(3/2))**(1/3)*besseli(-1/3, 2*z**(3/2)/3))/3'),
    ('airyaiprime(z).rewrite(hyper)', '3**(1/3)*z**2*hyper((), (5/3,), z**3/9)/(6*gamma(2/3)) - 3**(2/3)*hyper((), (1/3,), z**3/9)/(3*gamma(1/3))'),
    ('airybiprime(z).rewrite(besselj)', 'airybiprime(z)'),
    ('betainc(a, b, x, y).rewrite(hyper)', '(-x**a*hyper((a, 1 - b), (a + 1,), x) + y**a*hyper((a, 1 - b), (a + 1,), y))/a'),
    ('betainc_regularized(a, b, x, y).rewrite(Integral)', 'Integral(_t**(a - 1)*(1 - _t)**(b - 1), (_t, x, y))/Integral(_t**(a - 1)*(1 - _t)**(b - 1), (_t, 0, 1))'),
    ('zeta(s).rewrite(dirichlet_eta)', 'dirichlet_eta(s)/(1 - 2**(1 - s))'),
    ('zeta(s, x).rewrite(lerchphi)', 'lerchphi(1, s, x)'),
    ('polylog(s, z).rewrite(lerchphi)', 'z*lerchphi(z, s, 1)'),
    ('dirichlet_eta(s).rewrite(zeta)', 'Piecewise((log(2), Eq(s, 1)), ((1 - 2**(1 - s))*zeta(s), True))'),
    ('catalan(n).rewrite(binomial)', 'binomial(2*n, n)/(n + 1)'),
    ('catalan(n).rewrite(gamma)', '4**n*gamma(n + 1/2)/(sqrt(pi)*gamma(n + 2))'),
    ('catalan(n).rewrite(hyper)', 'hyper((-n, 1 - n), (2,), 1)'),
    ('catalan(n).rewrite(Product)', 'catalan(n)'),
    ('besselj(nu, z).rewrite(jn)', 'sqrt(2)*sqrt(z)*jn(nu - 1/2, z)/sqrt(pi)'),
    ('besselj(nu, z).rewrite(bessely)', 'besselj(nu, z)'),
    ('bessely(nu, z).rewrite(besselj)', 'bessely(nu, z)'),
    ('bessely(nu, z).rewrite(yn)', 'sqrt(2)*sqrt(z)*yn(nu - 1/2, z)/sqrt(pi)'),
    ('besseli(nu, z).rewrite(besselj)', 'exp(-I*pi*nu/2)*besselj(nu, z*exp_polar(I*pi/2))'),
    ('besselk(nu, z).rewrite(besseli)', 'besselk(nu, z)'),
    ('jn(nu, z).rewrite(besselj)', 'sqrt(2)*sqrt(pi)*sqrt(1/z)*besselj(nu + 1/2, z)/2'),
    ('jn(nu, z).rewrite(bessely)', '(-1)**nu*sqrt(2)*sqrt(pi)*sqrt(1/z)*bessely(-nu - 1/2, z)/2'),
    ('jn(nu, z).rewrite(yn)', '(-1)**nu*yn(-nu - 1, z)'),
    ('yn(nu, z).rewrite(besselj)', 'yn(nu, z)'),
    ('yn(2, z).rewrite(besselj)', '-sqrt(2)*sqrt(pi)*sqrt(1/z)*besselj(-5/2, z)/2'),
    ('yn(2, z).rewrite(bessely)', 'sqrt(2)*sqrt(pi)*sqrt(1/z)*bessely(5/2, z)/2'),
    ('yn(nu, z).rewrite(jn)', '(-1)**(nu + 1)*jn(-nu - 1, z)'),
    ('elliptic_k(m).rewrite(hyper)', 'pi*hyper((1/2, 1/2), (1,), m)/2'),
    ('elliptic_e(m).rewrite(hyper)', 'pi*hyper((-1/2, 1/2), (1,), m)/2'),
    ('elliptic_k(m).rewrite(meijerg)', 'meijerg(((1/2, 1/2), ()), ((0,), (0,)), -m)/2'),
    ('elliptic_f(z, m).rewrite(Integral)', 'Integral(1/sqrt(-m*sin(_t)**2 + 1), (_t, 0, z))'),
    ('elliptic_pi(n, z, m).rewrite(Integral)', 'Integral(1/(sqrt(-m*sin(_t)**2 + 1)*(-n*sin(_t)**2 + 1)), (_t, 0, z))'),
    ('elliptic_pi(1, z, m)', 'elliptic_f(z, m) + (sqrt(-m*sin(z)**2 + 1)*tan(z) - elliptic_e(z, m))/(1 - m)'),
    ('elliptic_pi(n, z, 0)', 'atanh(sqrt(n - 1)*tan(z))/sqrt(n - 1)'),
    ('elliptic_pi(n, z, n)', 'tan(z)/sqrt(-n*sin(z)**2 + 1) - (sqrt(-n*sin(z)**2 + 1)*tan(z) - elliptic_e(z, n))/(1 - n)'),
    ('stieltjes(-1)', 'zoo'),
    ('stieltjes(n)', 'stieltjes(n)'),
    ('stieltjes(0)', 'EulerGamma'),
    ('stieltjes(0, 2)', 'stieltjes(0, 2)'),
    ('stieltjes(S(1)/2)', 'zoo'),
    ('riemann_xi(oo)', 'nan'),
    ('riemann_xi(2)', 'pi/6'),
    ('hyper((3, 1, 2), [3, 4], x, evaluate=False)', 'hyper((1, 2, 3), (3, 4), x)'),
]


@pytest.mark.parametrize("source, expected", CASES)
def test_value(source, expected):
    assert str(eval(source, NS)) == expected


ERRORS = [
    ('Heaviside(x).fdiff(2)', ArgumentIndexError, 'Invalid operation with argument number 2 for Function Heaviside(x)'),
    ('DiracDelta(x).fdiff(2)', ArgumentIndexError, 'Invalid operation with argument number 2 for Function DiracDelta(x)'),
]


@pytest.mark.parametrize("source, kind, message", ERRORS)
def test_error(source, kind, message):
    with pytest.raises(kind) as info:
        eval(source, NS)
    assert str(info.value) == message


def test_singularity_derivative_form_matches_reference():
    # The multiple-derivative factor_terms(signsimp(...)) step gives the SymPy form.
    assert str(SingularityFunction(x, a, -2).rewrite(Heaviside)) == '-DiracDelta(a - x, 1)'


def test_documented_divergences_keep_correct_values():
    # SymPy returns 0 and nan; the order -1 polygamma diverges and
    # besselj(1, -I*y) = -I*besseli(1, y).
    assert polygamma(-1, oo) == oo
    assert besselj(1, -I*oo) == -oo*I


def test_stieltjes_numeric_single_argument():
    assert str(stieltjes(2).evalf(10)) == '-0.009690363193'
    assert str(stieltjes(1, 2).evalf(10)) == '-0.07281584548'
