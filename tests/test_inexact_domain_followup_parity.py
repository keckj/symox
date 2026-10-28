"""Inferred complex grounds, explicit conversions and inexact cancellation.

Compare accepted domains, coefficients and public return types with the local
SymPy 1.14 reference. The runtime implementation never imports SymPy.
"""
import pytest
import sympy as sp
import symox as sx


def namespace(module):
    ns = vars(module).copy()
    ns.update(zip(('x', 'y', 'z'), module.symbols('x y z')))
    return ns


def evaluate(module, source):
    return eval(source, namespace(module))


@pytest.mark.parametrize('source', [
    'Poly(x + 0.5*I)',
    'Poly(x + 0.1*y + 0.2*I)',
    'Poly(x*y + Float("0.2")*I*z + 2)',
    'Poly(sin(x) + Float("0.2")*I)',
    'Poly(x + Float("0.12345678901234567890123456789", 40)*I)',
    'Poly((x + Float("0.125")*I)**2)',
    # An explicitly requested I stays a polynomial generator over RR.
    'Poly(x + 0.5*I, x, I)',
    'Poly(x + 0.5*I, x)',
    'Poly(x + I/2)',
    'Poly(x + Float("0.2"))',
])
def test_inferred_complex_field_keeps_generators_and_coefficient_precision(source):
    actual, expected = evaluate(sx, source), evaluate(sp, source)
    assert tuple(map(str, actual.gens)) == tuple(map(str, expected.gens))
    assert str(actual.domain) == str(expected.domain)
    if hasattr(expected.domain, 'precision'):
        assert actual.domain.precision == expected.domain.precision
    assert [(m, sx.to_sympy(c)) for m, c in actual.terms()] == expected.terms()


@pytest.mark.parametrize('source', [
    'Poly(0.5*I)',
    'Poly(Float("0.1234567890123456789", 30)*I)',
])
def test_complex_constant_without_generators_still_needs_generators(source):
    for module in (sx, sp):
        with pytest.raises(module.polys.polyerrors.GeneratorsNeeded):
            evaluate(module, source)


@pytest.mark.parametrize('source', [
    'Poly(x**2 + 1, x, domain=CC).set_domain(QQ_I)',
    'Poly(x + Float("0.1")*I, x, domain=CC).set_domain(QQ_I)',
    'Poly(Float("0.125")*x*y + Float("0.2")*I, x, y, domain=CC).set_domain(QQ_I)',
    'Poly(Float("0.123456789012345678901234567890", 40)*x + I, x, domain=ComplexField(120)).set_domain(QQ_I)',
    'Poly(Float("1e-50", 40)*x + Float("0.1", 40)*I, x, domain=ComplexField(140)).set_domain(QQ_I)',
    'Poly(0, x, domain=CC).set_domain(QQ_I)',
    'Poly(x + Float("0.1", 40)*I).set_domain(QQ_I)',
])
def test_explicit_complex_to_gaussian_rational_conversion(source):
    actual, expected = evaluate(sx, source), evaluate(sp, source)
    assert str(actual.domain) == str(expected.domain) == 'QQ_I'
    assert [(m, sx.to_sympy(c)) for m, c in actual.terms()] == expected.terms()


@pytest.mark.parametrize('target', ['QQ', 'ZZ', 'ZZ_I'])
def test_other_exact_complex_conversions_keep_the_reference_rejection(target):
    source = 'Poly(x + 0.5*I, x, domain=CC).set_domain(%s)' % target
    for module in (sx, sp):
        with pytest.raises(module.polys.polyerrors.CoercionFailed):
            evaluate(module, source)


def test_complex_to_higher_precision_does_not_take_the_gaussian_rational_route():
    source = ('Poly(Float("0.123456789012345678901234567890", 40)*x + '
              'Float("0.987654321098765432109876543210", 40)*I, '
              'x, domain=ComplexField(120)).set_domain(ComplexField(160))')
    actual, expected = evaluate(sx, source), evaluate(sp, source)
    assert actual.domain.precision == expected.domain.precision == 160
    assert [(m, sx.to_sympy(c)) for m, c in actual.terms()] == expected.terms()


def test_real_symbolic_coefficients_are_evaluated_at_the_selected_precision():
    source = 'Poly(sqrt(2)*x + sin(1), x, domain=RealField(140))'
    actual, expected = evaluate(sx, source), evaluate(sp, source)
    assert actual.domain.precision == expected.domain.precision == 140
    assert [(m, sx.to_sympy(c)) for m, c in actual.terms()] == expected.terms()


CANCEL = [
    'cancel((1.5*x**2 - 1.5)/(x - 1), domain=RR)',
    'cancel((1.5*x**2 - 1.5)/(x - 1), domain="RR")',
    'cancel((x**2 + (I - 1)*x - I)/(x - 1), domain=CC)',
    'cancel((x**2 + (I - 1)*x - I)/(x - 1), domain="CC")',
    'cancel((1.5*y*x**2 - 1.5*y)/(x - 1), x, domain=RR[y])',
    'cancel((1.5*y*x**2 - 1.5*y)/(x - 1), x, domain=RR.frac_field(y))',
    'cancel((1.5*I*y*x**2 - 1.5*I*y)/(x - 1), x, domain=CC[y])',
    'cancel((1.5*I*y*x**2 - 1.5*I*y)/(x - 1), x, domain=CC.frac_field(y))',
    'cancel((1.5*x**2 + 1.5*x*z - 1.5*y*x - 1.5*y*z)/(x - y), x, z, domain=RR[y])',
    'cancel((Float("1.2345678901234567890123456789", 40)*x**2 - Float("1.2345678901234567890123456789", 40))/(x - 1), domain=RealField(140))',
    'cancel((1.5*x**2 - 1.5, x - 1), domain=RR)',
    'cancel((1.5*x**2 - 1.5, x - 1), domain=RR, polys=True)',
    'cancel((1.5*x**2 - 1.5, x - 1), domain=RR, polys=False)',
    'cancel((Poly(1.5*x**2 - 1.5, x), Poly(x - 1, x)), domain=RR)',
    'cancel((Poly(1.5*x**2 - 1.5, x), Poly(x - 1, x)), domain=RR, polys=False)',
    'cancel((x - 1)**2/(x - 1), domain=RR, expand=False)',
    'cancel(x/2 + 1, domain=RR)',
    'cancel((1, 2), domain=RR)',
    'cancel(Rational(1, 3), domain=RR)',
    'cancel(Eq(x, 1), domain=RR)',
    'cancel(Poly(x + 1, x), domain=RR)',
    'cancel((1.5*x**2 - 1.5)/(1 - x), domain=RR, _signsimp=False)',
    'cancel(Piecewise(((x**2 - 1)/(x - 1), x > 0), (x, True)), domain=RR)',
    'cancel(Symbol("A", commutative=False)*(x**2 - 1)/(x - 1), domain=RR)',
    'cancel((Poly(x**2 - 1, x), Poly(x - 1, x)), domain=RR)',
    'cancel((Poly(x**2 - 1, x, domain=QQ), Poly(x - 1, x)), domain=CC)',
    'cancel((x + I, x + 1), domain=CC, polys=True)',
    'cancel((x, x + 1), domain=RR, gens=(x,))',
    'cancel((x**2 - 1)/(x - 1), domain="R")',
    'cancel((x**2 - 1)/(x - 1), domain="RR_100")',
    'cancel((x**2 - 1)/(x - 1), domain="CC_100")',
    'cancel(sin(1), domain=RR)',
    'cancel(sqrt(2), domain=RR)',
    'cancel((sin(1), 2), domain=RR)',
    'cancel(sqrt(2)*I, domain=RR)',
    'cancel(I, domain=CC)',
    'cancel((Piecewise(((x**2 - 1)/(x - 1), x > 0), (x, True)), x + 1), domain=RR)',
    'cancel((Symbol("A", commutative=False)*x + Symbol("A", commutative=False), x + 1), domain=RR)',
]


def converted(value):
    if isinstance(value, sx.Poly):
        return ('Poly', tuple(map(str, value.gens)), str(value.domain),
                [(m, sx.to_sympy(c)) for m, c in value.terms()])
    if isinstance(value, sp.Poly):
        return ('Poly', tuple(map(str, value.gens)), str(value.domain), value.terms())
    if isinstance(value, tuple):
        return tuple(converted(v) for v in value)
    return sx.to_sympy(value) if isinstance(value, sx.Basic) else value


@pytest.mark.parametrize('source', CANCEL)
def test_cancel_explicit_inexact_domains_and_public_return_types(source):
    assert converted(evaluate(sx, source)) == converted(evaluate(sp, source))


@pytest.mark.parametrize('source, error', [
    ('cancel((x,), domain=RR)', 'ValueError'),
    ('cancel((x, x, x), domain=RR)', 'ValueError'),
    ('cancel(x, domain=RR, frac=True)', 'FlagError'),
    ('cancel(x, domain=RR, auto=True)', 'FlagError'),
    ('cancel(x, domain=RR, field=True)', 'OptionError'),
    ('cancel(x, domain=RR, greedy=False)', 'OptionError'),
    ('cancel((Poly(x**2 - 1, x), Poly(x - 1, x)), x, domain=RR)', 'OptionError'),
    ('cancel((x**2 - 1, x - 1), x, domain=RR, polys=True)', 'OptionError'),
])
def test_cancel_input_and_option_errors(source, error):
    for module in (sx, sp):
        exception = (ValueError if error == 'ValueError'
                     else getattr(module.polys.polyerrors, error))
        with pytest.raises(exception):
            evaluate(module, source)


@pytest.mark.parametrize('domain', ['CC[y]', 'CC.frac_field(y)'])
@pytest.mark.parametrize('method', ['is_positive', 'is_negative', 'is_nonpositive', 'is_nonnegative'])
def test_complex_composite_ground_has_no_order(domain, method):
    source = '(lambda K: K.%s(K.convert(1.5*I*y)))(%s)' % (method, domain)
    assert evaluate(sx, source) is evaluate(sp, source) is False


def test_high_precision_complex_cancel_does_not_fabricate_coefficient_digits():
    # SymPy reconstructs each mpc part through the default 53-bit context
    # during the exact-ground gcd. It loses digits despite returning 140-bit
    # Floats. Arithmetic must retain the selected field's precision instead.
    x = sx.Symbol('x')
    a = sx.Float('1.2345678901234567890123456789', 40)
    numerator = a*sx.I*x**2 - a*sx.I
    actual = sx.cancel(numerator/(x - 1), domain=sx.ComplexField(140))
    expected = sx.to_sympy(a)*(sp.I*sp.Symbol('x') + sp.I)
    assert sp.expand(sx.to_sympy(actual) - expected) == 0
    F, G = sx.Poly(numerator, x, domain=sx.ComplexField(140)), sx.Poly(x-1, x, domain=sx.ComplexField(140))
    common, p, q = F.cofactors(G)
    # Float object equality also compares _prec. Prove equality of the
    # stored binary values while requiring the selected 140-bit output.
    leading = sx.to_sympy(p.LC())
    assert sp.Rational(sp.re(leading)) == 0
    assert sp.Rational(sp.im(leading)) == sp.Rational(sx.to_sympy(a))
    assert p.domain.precision == 140
    assert {f._prec for f in leading.atoms(sp.Float)} == {140}
    assert sx.to_sympy(common.as_expr()*p.as_expr() - F.as_expr()).expand() == 0
    assert sx.to_sympy(common.as_expr()*q.as_expr() - G.as_expr()).expand() == 0
