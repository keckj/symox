"""Mellin gamma moments retain their exact convergence conditions."""
import pytest
import subprocess
import sys
import symox as sx
import sympy as sp


def same_value(actual, expected, rate):
    actual = actual.to_sympy()
    if actual == expected:
        return
    if rate == '1+I':
        # The reference keeps a polar lift of 1+I. Its principal logarithm
        # is exactly log(2)/2+I*pi/4; there is no winding in this right-half-
        # plane constant. Compare that exact branch in exponential form.
        actual = actual.replace(lambda v: v.is_Pow and v.base == 1+sp.I,
            lambda v: sp.exp(v.exp*(sp.log(2)/2+sp.I*sp.pi/4)))
        expected = sp.unpolarify(sp.powdenest(expected, polar=True))
        assert sp.simplify((actual-expected).rewrite(sp.exp)) == 0
    else:
        assert sp.simplify(actual-expected) == 0


@pytest.mark.parametrize('degree', [1, 2, 4])
@pytest.mark.parametrize('rate', ['1', '2', 'positive', '1+I'])
@pytest.mark.parametrize('conds', ['piecewise', 'none', 'separate'])
@pytest.mark.parametrize('method', [False, True])
def test_symbolic_moments_and_real_part_conditions(degree, rate, conds, method):
    x, a = sx.symbols('x a')
    beta = sx.Symbol('b', positive=True) if rate == 'positive' else sx.sympify(rate)
    expression = x**a*sx.exp(-beta*x**degree)
    function = expression.integrate if method else lambda *limits, **kwargs: sx.integrate(expression, *limits, **kwargs)
    actual = function((x, 0, sx.oo), conds=conds)
    expected = sp.integrate(expression.to_sympy(), (x.to_sympy(), 0, sp.oo), conds=conds)
    if conds == 'separate':
        assert isinstance(actual, tuple) and len(actual) == 2
        value, condition = actual
        assert condition.to_sympy() == expected[1]
        same_value(value, expected[0], rate)
    elif conds == 'piecewise':
        converted = actual.to_sympy()
        assert isinstance(converted, sp.Piecewise)
        assert converted.args[0][1] == expected.args[0][1]
        same_value(actual.args[0][0], expected.args[0][0], rate)
        assert converted.args[1] == expected.args[1]
    else:
        same_value(actual, expected, rate)


@pytest.mark.parametrize('exponent', ['I', '-1/2+I', '3'])
def test_provably_convergent_complex_powers_need_no_piecewise(exponent):
    x = sx.Symbol('x');p = sx.sympify(exponent)
    expression = x**p*sx.exp(-x)
    actual = sx.integrate(expression, (x, 0, sx.oo))
    assert not actual.is_Integral and not actual.is_Piecewise
    expected = sp.gamma(p.to_sympy()+1)
    assert sp.simplify(actual.to_sympy()-expected) == 0


@pytest.mark.parametrize('conds', ['piecewise', 'none', 'separate'])
def test_divergent_complex_moment_does_not_become_analytic_continuation(conds):
    x = sx.Symbol('x')
    expression = x**(-2+sx.I)*sx.exp(-x)
    actual = sx.integrate(expression, (x, 0, sx.oo), conds=conds)
    assert actual.is_Integral


def test_integral_doit_uses_the_same_default_conditions():
    x, a = sx.symbols('x a')
    expression = x**a*sx.exp(-x)
    assert sx.Integral(expression, (x, 0, sx.oo)).doit() == sx.integrate(expression, (x, 0, sx.oo))


def test_finite_and_indefinite_integrals_keep_scalar_results_with_separate():
    x = sx.Symbol('x')
    assert sx.integrate(x*x, (x, 0, 1), conds='separate') == sx.Rational(1, 3)
    assert sx.integrate(x*x, x, conds='separate') == x**3/3
    assert sx.integrate(sx.exp(-x), (x, 0, sx.oo), conds='separate') == (sx.Integer(1), sx.true)


@pytest.mark.parametrize('conds', ['invalid', None, True])
def test_invalid_condition_modes_are_rejected(conds):
    x = sx.Symbol('x')
    with pytest.raises(ValueError):
        sx.integrate(x, x, conds=conds)


@pytest.mark.parametrize('container', [tuple, list, sx.Tuple])
def test_limits_accept_python_and_native_sequences(container):
    x, a = sx.symbols('x a')
    limit = container(x, 0, sx.oo) if container is sx.Tuple else container([x, 0, sx.oo])
    assert sx.integrate(x**a*sx.exp(-x), limit, conds='separate') == (sx.gamma(a+1), sx.re(a)>-1)


@pytest.mark.parametrize('conds', ['none', 'separate'])
def test_reversed_improper_limits_keep_the_orientation(conds):
    x, a = sx.symbols('x a')
    actual = sx.integrate(x**a*sx.exp(-x), (x, sx.oo, 0), conds=conds)
    expected = (-sx.gamma(a+1), sx.re(a)>-1) if conds == 'separate' else -sx.gamma(a+1)
    assert actual == expected


def test_shifted_exponential_keeps_its_free_coefficient():
    x, a, c = sx.symbols('x a c')
    actual = sx.integrate(x**a*sx.exp(c-x), (x, 0, sx.oo), conds='separate')
    assert actual == (sx.exp(c)*sx.gamma(a+1), sx.re(a)>-1)


def test_cancelling_divergent_terms_keep_one_original_integral_fallback():
    x, a = sx.symbols('x a')
    expression = (sx.exp(-x)-sx.exp(-2*x))*x**a
    result = sx.integrate(expression, (x, 0, sx.oo))
    assert result.is_Piecewise
    # This boundary integral converges by cancellation. The elementary gamma
    # rule alone cannot decide it and must retain the original convergent
    # integral, rather than produce gamma(0)-gamma(0) or infinity-infinity.
    boundary = result.subs(a, -1)
    assert boundary.is_Integral
    assert boundary.function == (sx.exp(-x)-sx.exp(-2*x))/x


def test_separate_is_rejected_for_multiple_limits_like_reference():
    x, y = sx.symbols('x y')
    with pytest.raises(ValueError, match='multiple integrals'):
        sx.integrate(sx.exp(-x-y), (x, 0, sx.oo), (y, 0, sx.oo), conds='separate')


def test_large_power_moments_are_sparse_and_do_not_import_reference():
    script = '''
import sys
class NoReference:
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'sympy' or fullname.startswith('sympy.'):
            raise ImportError('reference forbidden')
sys.meta_path.insert(0, NoReference())
import symox as sx
x,a=sx.symbols('x a')
q=10**6
value,condition=sx.integrate(x**a*sx.exp(3-(1+sx.I)*x**q),[x,0,sx.oo],conds='separate')
s=(a+1)/q
assert value==(1-sx.exp_polar(-sx.I*sx.pi/2))**(-s)*sx.exp(3)*sx.gamma(s)/q
assert condition==(sx.re(a)>-1)
assert sx.integrate(x**a*sx.exp(-x**q),[x,0,sx.oo]).is_Piecewise
'''
    subprocess.run([sys.executable, '-c', script], check=True, timeout=20)
