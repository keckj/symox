"""Precision of expm1/log1p and mpmath lambdify (wave-4 replay, category e).

Accuracy is measured against mpmath at 400 digits, not against SymPy, which
loses every digit in several of these cases.
"""
import mpmath
import pytest

import symox as S
from symox.codegen.cfunctions import expm1, log1p
from symox.printing.pycode import MpmathPrinter

x, y = S.symbols('x y')


def relative_error(value, reference):
    with mpmath.workdps(400):
        return abs(mpmath.mpf(str(value)) - reference) / abs(reference)


@pytest.mark.parametrize('dps', [15, 30, 60, 120])
@pytest.mark.parametrize('k', [10, 50, 99, 200])
def test_expm1_and_log1p_have_no_cancellation(dps, k):
    t = S.Rational(1, 10**k)
    with mpmath.workdps(400):
        small = mpmath.mpf(1) / mpmath.mpf(10)**k
        e_ref, l_ref = mpmath.expm1(small), mpmath.log1p(small)
    bound = 10.0**(1 - dps)
    for value, reference in [(expm1(t), e_ref), (log1p(t), l_ref), (S.expand_log(log1p(t)), l_ref),
                             (expm1(S.Float(t, dps + 10)), e_ref), (log1p(S.Float(t, dps + 10)), l_ref),
                             (S.expand_log(log1p(S.Float(t, dps + 10))), l_ref),
                             (expm1(x).subs(x, t), e_ref), (log1p(x).subs(x, t), l_ref)]:
        assert relative_error(S.N(value, dps), reference) < bound


def test_reference_examples():
    assert '%.0e' % expm1(1e-99).evalf() == '1e-99'
    assert '%.0e' % S.expand_log(log1p(1e-99)).evalf() == '1e-99'
    # SymPy evaluates log(1 + 1e-99) to 0; the exact form is kept here.
    assert '%.0e' % log1p(1e-99).evalf() == '1e-99'
    assert log1p(1e-99) == S.log(S.Rational(S.Float(1e-99)) + 1)


def test_expm1_log1p_symbolic_behaviour():
    assert expm1(x).diff(x) == S.exp(x)
    assert log1p(x).diff(x) == 1/(x + 1)
    assert expm1(0) == 0 and expm1(S.log(3)) == 2
    assert log1p(2) == S.log(3) and log1p(0) == 0
    assert str(expm1(2)) == 'expm1(2)'
    assert S.N(expm1(2), 30) == S.N(S.exp(2) - 1, 30)
    assert S.expand_func(expm1(x)) == S.exp(x) - 1
    assert expm1(x).rewrite(S.exp) == S.exp(x) - 1
    assert log1p(x).rewrite(S.log) == S.log(x + 1)
    assert expm1(S.Symbol('r', real=True)).is_real is True
    assert log1p(S.Symbol('p', positive=True)).is_positive is True
    assert isinstance(expm1(x), expm1)
    assert S.ccode(expm1(x) + log1p(x)) == 'expm1(x) + log1p(x)'
    assert S.lambdify(x, expm1(x) + log1p(x))(1e-30) == pytest.approx(2e-30, rel=1e-15)
    z = expm1(S.Float('1e-40', 30) * (1 + S.I))
    assert abs(complex(z) - complex(1e-40, 1e-40)) < 1e-55


def test_cancellation_sensitive_forms_are_at_least_reference_accurate():
    t = S.Rational(1, 10**99)
    with mpmath.workdps(400):
        small = mpmath.mpf(1) / mpmath.mpf(10)**99
        cases = [(S.cos(t) - 1, mpmath.cos(small) - 1), (S.cosh(t) - 1, mpmath.cosh(small) - 1),
                 (S.log(1 + t) - t, mpmath.log1p(small) - small), (S.sinc(t) - 1, mpmath.sin(small)/small - 1),
                 (S.tan(t) - t, mpmath.tan(small) - small), (S.exp(t) - 1 - t, mpmath.expm1(small) - small),
                 (S.Si(t), mpmath.si(small)), (S.erf(t), mpmath.erf(small)), (S.sqrt(1 + t) - 1, mpmath.sqrt(1 + small) - 1)]
    for expr, reference in cases:
        for dps in (15, 30, 60):
            assert relative_error(S.N(expr, dps), reference) < 10.0**(2 - dps)


def test_lambdify_mpmath_keeps_the_working_precision():
    with mpmath.workdps(40):
        value = S.lambdify(x, S.sin(x), 'mpmath')(1)
        assert isinstance(value, mpmath.mpf)
        assert mpmath.nstr(value, 40) == mpmath.nstr(mpmath.sin(1), 40)
        f = S.lambdify(x, S.pi*S.E + S.sqrt(x) + S.Rational(1, 3) + S.Float('1.23456789012345678901234567890', 30)*x,
                       modules=['mpmath'])
        expected = mpmath.pi*mpmath.e + mpmath.sqrt(2) + mpmath.mpf(1)/3 + mpmath.mpf('1.23456789012345678901234567890')*2
        assert abs(f(2) - expected) < mpmath.mpf(10)**-29
        g = S.lambdify(x, x**S.Rational(3, 2) + S.EulerGamma + S.LambertW(x) + S.Si(x), 'mpmath')
        assert abs(g(2) - (mpmath.mpf(2)**1.5 + mpmath.euler + mpmath.lambertw(2) + mpmath.si(2))) < mpmath.mpf(10)**-38
        assert S.lambdify(x, expm1(x) + log1p(x), 'mpmath')(mpmath.mpf('1e-30')) == pytest.approx(2e-30, rel=1e-30)
        # Namespace precedence, structured arguments and outputs.
        assert S.lambdify(x, S.sin(x), [{'sin': mpmath.cos}, 'mpmath'])(0) == 1
        assert S.lambdify([(x, y)], x + y, 'mpmath')((1, 2)) == 3
        assert S.lambdify(x, [x, S.pi], 'mpmath')(1)[1] == mpmath.pi
        assert isinstance(S.lambdify(x, S.Matrix([[x, 1]]), 'mpmath')(2), mpmath.matrix)
        assert S.lambdify(x, S.Piecewise((x, x > 0), (0, True)), 'mpmath')(-1) == 0
        # The math namespace comes first: SymPy also computes doubles there.
        assert isinstance(S.lambdify(x, S.sin(x), ['math', 'mpmath'])(1), float)


def test_mpmath_printer_matches_reference():
    printer = MpmathPrinter()
    assert printer.doprint(S.Float('1.23456789012345678901234567890', 30)*x) == \
        'mpmath.mpf((0, 6260002907749995053341990378197, -102, 103))*x'
    assert printer.doprint(S.sqrt(x) + S.Rational(1, 3)) == 'mpmath.sqrt(x) + mpmath.mpf(1)/mpmath.mpf(3)'
    assert printer.doprint(x**S.Rational(3, 2)) == 'x**(mpmath.mpf(3)/mpmath.mpf(2))'
    assert printer.doprint(1/S.sqrt(x)) == '1/mpmath.sqrt(x)'
    assert printer.doprint(S.oo) == 'mpmath.inf'
    assert printer.doprint(S.EulerGamma) == 'mpmath.euler'
