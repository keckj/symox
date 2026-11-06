"""The native ``RootSum`` constructor (wave 43, lot ratint2).

``RootSum.__new__`` after the argument checks, ``_transform``
(``preprocess_roots``), ``new``, ``_rational_case``, ``doit`` and
``_eval_evalf`` compute in Rust (``crates/symox-py/src/root_sum.rs`` on
``symox::poly::rootsum`` and the closed forms of ``symox::poly::roots``).
The expected strings are SymPy 1.14 outputs, hard-coded. The tests that
replace the Python helpers by a failing stub prove that the native path is
taken.
"""
import pytest

import symox as sx
from symox import Lambda, RootSum, log, exp, sqrt, Rational
from symox import _core, _public_roots
from symox.polys import polyroots, polyfuncs
from symox.polys.polyerrors import MultivariatePolynomialError

x, y, z, a = sx.symbols('x y z a')

# (description, constructor, SymPy 1.14 str)
REFERENCE = [
    ('identity', lambda: RootSum(x**3 + x + 1), '0'),
    ('1/(x+1)', lambda: RootSum(x**3 + x + 1, Lambda(x, 1/(x + 1))), '4'),
    ('x**2', lambda: RootSum(x**3 + x + 1, Lambda(x, x**2)), '-2'),
    ('x/(x**2+2)', lambda: RootSum(x**3 + x + 1, Lambda(x, x/(x**2 + 2))), '5/3'),
    ('log', lambda: RootSum(x**3 + x + 1, Lambda(x, log(x))), 'RootSum(x**3 + x + 1, Lambda(x, log(x)))'),
    ('2+3*y*log', lambda: RootSum(x**3 + x + 1, Lambda(y, 2 + 3*y*log(z - y))),
     '3*RootSum(x**3 + x + 1, Lambda(y, y*log(-y + z))) + 6'),
    ('exp', lambda: RootSum(x**3 + x + 1, Lambda(x, exp(x))), 'RootSum(x**3 + x + 1, Lambda(x, exp(x)))'),
    ('Function f', lambda: RootSum(x**3 + x + 1, sx.Function('f')), 'RootSum(x**3 + x + 1, Lambda(x, f(x)))'),
    ('constant', lambda: RootSum(x**3 + x + 1, Lambda(x, 7)), '21'),
    ('quintic x**3', lambda: RootSum(x**5 - x + 1, Lambda(x, x**3)), '0'),
    ('quintic 1/x', lambda: RootSum(x**5 - x + 1, Lambda(x, 1/x)), '1'),
    ('product', lambda: RootSum((x**2 + 1)*(x**3 + x + 1), Lambda(x, log(x))),
     'RootSum(x**2 + 1, Lambda(x, log(x))) + RootSum(x**3 + x + 1, Lambda(x, log(x)))'),
    ('quadratic=True', lambda: RootSum(x**2 + 3*x + 1, Lambda(x, log(x)), quadratic=True),
     'log(3/2 - sqrt(5)/2) + log(sqrt(5)/2 + 3/2) + 2*I*pi'),
    ('linear factor', lambda: RootSum((x - 2)*(x**3 + x + 1), Lambda(x, log(x))),
     'RootSum(x**3 + x + 1, Lambda(x, log(x))) + log(2)'),
    ('repeated factor', lambda: RootSum((x**3 + x + 1)**2, Lambda(x, log(x))),
     '2*RootSum(x**3 + x + 1, Lambda(x, log(x)))'),
    ('auto=False', lambda: RootSum(x**3 + x + 1, Lambda(x, 1/(x + 1)), auto=False),
     'RootSum(x**3 + x + 1, Lambda(x, 1/(x + 1)))'),
    ('scaled', lambda: RootSum(4*x**3 + 2*x + 3, Lambda(x, log(x))), 'RootSum(4*x**3 + 2*x + 3, Lambda(x, log(x)))'),
    ('integer basis', lambda: RootSum(x**3 + 8*x + 64, Lambda(x, log(x))), 'RootSum(x**3 + 2*x + 8, Lambda(x, log(2*x)))'),
    ('integer basis rational', lambda: RootSum(x**3 + 8*x + 64, Lambda(x, 1/(x + 1))), '-1/5'),
    ('sqrt2 log', lambda: RootSum(x**3 + sqrt(2)*x + 1, Lambda(x, log(x))),
     'RootSum(x**3 + sqrt(2)*x + 1, Lambda(x, log(x)))'),
    ('sqrt2 rational', lambda: RootSum(x**3 + sqrt(2)*x + 1, Lambda(x, x**2)), '-2*sqrt(2)'),
    ('fraction coefficients', lambda: RootSum(x**3 + x/2 + Rational(1, 3), Lambda(x, log(x))),
     'RootSum(6*x**3 + 3*x + 2, Lambda(x, log(x)))'),
    ('parametric gen', lambda: RootSum(x**3 + a*x + 1, Lambda(x, a/(x + 1)), x), 'a + 3'),
    ('parametric gen rational', lambda: RootSum(x**3 + a*x + 1, Lambda(x, a*x**2 + x/(x + a)), x),
     '-(2*a**5 + 2*a**4 - 4*a**2 + 3)/(a**3 + a**2 - 1)'),
    ('sqrt2 summand', lambda: RootSum(x**3 + x + 1, Lambda(x, sqrt(2)*x/(x + sqrt(2)))),
     'sqrt(2)*(-3 + 2*sqrt(2))/(-1 + 3*sqrt(2))'),
    ('doit quartic', lambda: RootSum(x**4 + 1, Lambda(x, log(x))).doit(), '0'),
    ('doit quintic stays', lambda: RootSum(x**5 - x + 1, Lambda(x, log(x))).doit(),
     'RootSum(x**5 - x + 1, Lambda(x, log(x)))'),
    ('evalf cubic', lambda: RootSum(x**3 + x + 1, Lambda(x, log(x - 1))).evalf(),
     '1.09861228866811 + 3.14159265358979*I'),
    ('evalf quintic', lambda: RootSum(x**5 - x + 1, Lambda(x, log(x + 2))).evalf(30),
     '3.36729582998647402718327203236'),
    ('diff', lambda: sx.diff(RootSum(x**3 + x + 1, Lambda(y, y*log(z - y))), z), '-(2*z + 3)/(z**3 + z + 1)'),
    ('subs', lambda: RootSum(x**3 + x + 1, Lambda(y, log(z - y))).subs(z, 2),
     'RootSum(x**3 + x + 1, Lambda(y, log(2 - y)))'),
]


@pytest.mark.parametrize('name,build,expected', REFERENCE, ids=[c[0] for c in REFERENCE])
def test_root_sum_matches_reference(name, build, expected):
    assert str(build()) == expected


@pytest.mark.parametrize('build', [
    lambda: RootSum(x**3 + a*x + a**3, Lambda(x, log(x))),
    lambda: RootSum(x**4 + a*x + 1, Lambda(x, 1/(x - a))),
    lambda: RootSum(x**3 + y*x + 1 + z*y),
])
def test_root_sum_multivariate_error(build):
    with pytest.raises(MultivariatePolynomialError, match='only univariate polynomials are allowed'):
        build()


def test_root_sum_bad_function_error():
    with pytest.raises(ValueError, match='expected a univariate function'):
        RootSum(x**3 + x + 1, Lambda((x, y), x + y))


def test_raw_node_records_polynomial_and_auto():
    r = RootSum(x**3 + x + 1, Lambda(x, 1/(x + 1)), auto=False)
    assert r.auto is False
    assert str(r.poly) == "PurePoly(x**3 + x + 1, x, domain='ZZ')"
    assert r.args == (x**3 + x + 1, Lambda(x, 1/(x + 1)), x)


def _fail(*args, **kwargs):
    raise AssertionError('the Python fallback ran')


def test_constructor_native_path(monkeypatch):
    monkeypatch.setattr(polyroots, 'preprocess_roots', _fail)
    monkeypatch.setattr(polyroots, 'roots_linear', _fail)
    monkeypatch.setattr(polyroots, 'roots_quadratic', _fail)
    monkeypatch.setattr(polyfuncs, 'symmetrize', _fail)
    monkeypatch.setattr(polyfuncs, 'viete', _fail)
    monkeypatch.setattr(_public_roots.RootSum, '_rational_case', classmethod(_fail))
    monkeypatch.setattr(_public_roots, '_pure_factors', _fail)
    assert str(RootSum(x**3 + 8*x + 64, Lambda(x, 1/(x + 1)))) == '-1/5'
    assert str(RootSum(x**3 + a*x + 1, Lambda(x, a*x**2 + x/(x + a)), x)) == \
        '-(2*a**5 + 2*a**4 - 4*a**2 + 3)/(a**3 + a**2 - 1)'
    assert str(RootSum((x - 2)*(x**2 + 3*x + 1), Lambda(x, log(x)), quadratic=True)) == \
        str(sx.log(2) + RootSum(x**2 + 3*x + 1, Lambda(x, log(x)), quadratic=True))
    assert str(RootSum(x**3 + x + 1, Lambda(x, log(x)))) == 'RootSum(x**3 + x + 1, Lambda(x, log(x)))'
    # `RootSum.new` (the derivative) is native too.
    assert str(sx.diff(RootSum(x**3 + x + 1, Lambda(y, y*log(z - y))), z)) == '-(2*z + 3)/(z**3 + z + 1)'


def test_doit_and_evalf_native_path(monkeypatch):
    r = RootSum(x**4 + 1, Lambda(x, log(x)))
    e = RootSum(x**3 + x + 1, Lambda(x, log(x - 1)))
    monkeypatch.setattr(polyroots, 'roots', _fail)
    monkeypatch.setattr(sx.PurePoly, 'nroots', _fail)
    assert str(r.doit()) == '0'
    assert str(e.evalf()) == '1.09861228866811 + 3.14159265358979*I'


def test_native_functions_decline_without_poly():
    assert _core._root_sum_transform(x) is None
    assert _core._root_sum_build(1, x, x, x) is None
