"""Exact rational-reduction certificates, domains and complete fallbacks."""
import pytest
import symox as S
import sympy as R
from symox import _core
from symox.integrals.rationaltools import ratint_ratpart
from sympy.integrals.rationaltools import ratint_ratpart as reference_reduction

x = S.Symbol('x')
rx = R.Symbol('x')


@pytest.mark.parametrize('variant', range(24))
def test_exact_repeated_poles_match_reference_and_polynomial_certificate(variant):
    # Distinct rational coefficients, varying multiplicity and square-free part.
    a, b, multiplicity = variant % 5 + 1, variant % 7 + 2, variant % 4 + 1
    denominator = (x - S.Rational(a, b))**multiplicity * (x*x + x + b)
    numerator = (variant + 1)*x*x + S.Rational(2*variant - 3, b)*x + 1
    f, g = S.Poly(numerator, x, domain='QQ'), S.Poly(denominator, x, domain='QQ')
    result = _core._ratint_ratpart(f, g, x)
    assert result is not None
    expected = reference_reduction(R.Poly(str(numerator), rx), R.Poly(str(denominator), rx), rx)
    actual = tuple(R.sympify(str(value), locals={'x': rx}) for value in result)
    assert all(R.cancel(a - b) == 0 for a, b in zip(actual, expected))
    assert R.cancel(R.diff(actual[0], rx) + actual[1] - R.sympify(str(numerator/denominator))) == 0
    log_denominator = R.Poly(R.fraction(R.cancel(actual[1]))[1], rx)
    assert R.gcd(log_denominator, log_denominator.diff()).degree() == 0
    assert ratint_ratpart(f, g, x) == result


@pytest.mark.parametrize('domain', ['ZZ', 'QQ'])
def test_public_output_and_square_free_case(domain):
    f, g = S.Poly(36, x, domain=domain), S.Poly(x**5 - 2*x**4 - 2*x**3 + 4*x*x + x - 2, x, domain=domain)
    rational, logarithmic = ratint_ratpart(f, g, x)
    assert str(rational) == '(12*x + 6)/(x**2 - 1)'
    assert str(logarithmic) == '12/(x**2 - x - 2)'
    assert ratint_ratpart(S.Poly(1, x), S.Poly(2*x + 1, x), x) == (0, 1/(2*x + 1))
    assert ratint_ratpart(S.Poly(0, x), S.Poly(x*x, x), x) == (0, 0)


@pytest.mark.parametrize('domain', ['EX', 'QQ(a)', 'RR'])
def test_other_domains_decline_and_preserve_complete_host_reduction(domain):
    a = S.Symbol('a')
    numerator = 1 if domain == 'RR' else a
    f, g = S.Poly(numerator, x, domain=domain), S.Poly((x + 1)**2, x, domain=domain)
    assert _core._ratint_ratpart(f, g, x) is None
    rational, logarithmic = ratint_ratpart(f, g, x)
    assert S.cancel(rational.diff(x) + logarithmic - f.as_expr()/g.as_expr()) == 0


def test_route_boundaries_are_refusals_without_densifying_huge_coordinates():
    one = S.Poly(1, x)
    assert _core._ratint_ratpart(one, S.Poly(x**129 + 1, x), x) is None
    assert _core._ratint_ratpart(one, S.Poly(x**(2**31) + 1, x), x) is None
    assert _core._ratint_ratpart(one, one, x) is None
    assert _core._ratint_ratpart(one, S.Poly(0, x), x) is None
    assert _core._ratint_ratpart(S.Poly(x*x, x), S.Poly(x + 1, x), x) is None
    assert _core._ratint_ratpart(one, S.Poly(x*x, x), S.Symbol('other')) is None


def test_custom_polynomial_subclass_keeps_its_protocol():
    class CustomPoly(S.Poly):
        pass
    assert _core._ratint_ratpart(CustomPoly(1, x), S.Poly(x*x, x), x) is None


def test_dummy_generator_and_large_exact_coefficients():
    d = S.Dummy('z')
    f, g = S.Poly(2**257 + 1, d), S.Poly((3*d + 2)**3, d)
    result = _core._ratint_ratpart(f, g, d)
    assert result is not None
    assert S.cancel(result[0].diff(d) + result[1] - f.as_expr()/g.as_expr()) == 0


def test_general_expression_inputs_still_use_public_polynomial_conversion():
    assert ratint_ratpart(1, x*x, x) == (-1/x, 0)


def test_method_changes_decline_before_bypassing_custom_protocol(monkeypatch):
    f, g = S.Poly(1, x), S.Poly(x*x, x)
    original = S.Poly.diff
    try:
        monkeypatch.setattr(S.Poly, 'diff', lambda *args, **kwargs: original(*args, **kwargs))
    except TypeError:
        pytest.skip('the native polynomial class does not allow changing methods')
    assert _core._ratint_ratpart(f, g, x) is None
