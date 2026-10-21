"""Improper definite integrals never report a value without a convergence proof.

An antiderivative evaluated at an infinite bound gives a value only when the
limit is settled: a signed infinity (a proven divergence) or a finite value
free of every infinity. Otherwise the integral is returned unevaluated, as
the reference does when it cannot prove convergence.
"""
import pytest
import symox as sx


x, t = sx.symbols('x t')
s, a = sx.symbols('s a')
p = sx.Symbol('p', positive=True)


@pytest.mark.parametrize('integrand, bounds', [
    (lambda: sx.exp(-s*t), (0, sx.oo)),
    (lambda: t*sx.exp(-s*t), (0, sx.oo)),
    (lambda: sx.sin(x), (0, sx.oo)),
    (lambda: sx.cos(s*x), (0, sx.oo)),
    (lambda: x**s, (1, sx.oo)),
])
def test_unproven_limit_keeps_the_integral(integrand, bounds):
    f = integrand()
    var = t if f.has(t) else x
    result = sx.integrate(f, (var, *bounds))
    if result.has(sx.Integral):
        return
    # A value is acceptable only together with its convergence condition.
    assert isinstance(result, sx.Piecewise)
    assert result.args[-1][0].has(sx.Integral)


def test_settled_limits_still_evaluate():
    assert sx.integrate(sx.exp(-p*t), (t, 0, sx.oo)) == 1/p
    assert sx.integrate(1/(1+x**2), (x, -sx.oo, sx.oo)) == sx.pi
    assert sx.integrate(1/x, (x, 1, sx.oo)) == sx.oo
    assert sx.integrate(sx.exp(x), (x, 0, sx.oo)) == sx.oo


def test_hyperbolic_laplace_kernel_keeps_the_meijer_result_and_condition():
    # The reference keeps the Meijer G evaluation, whose Piecewise has a
    # G-function in its last piece, rather than the exp rewrite.
    result = sx.integrate(sx.exp(-s*t)*sx.cosh(t), (t, 0, sx.oo))
    assert isinstance(result, sx.Piecewise)
    (value, cond), (fallback, otherwise) = result.args
    assert cond == (s**2 > 1)
    assert fallback == sx.Integral(sx.exp(-s*t)*sx.cosh(t), (t, 0, sx.oo))
    assert otherwise == sx.true
    assert value.has(sx.meijerg)
    assert sx.simplify(value.subs(s, 2)) == sx.Rational(2, 3)
