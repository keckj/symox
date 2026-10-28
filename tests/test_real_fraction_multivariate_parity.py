"""Complete real fraction-field gcd/cofactor/cancel routes and conventions."""
import pytest
import sympy as sp
import symox as sx
from mpmath.libmp import dps_to_prec, prec_to_dps


def evaluate(module, source):
    namespace = vars(module).copy()
    namespace.update(zip(('x', 'z', 'w', 'y', 't'), module.symbols('x z w y t')))
    return eval(source, namespace)


def exact(value):
    if isinstance(value, sx.Basic):
        value = sx.to_sympy(value)
    value = sp.sympify(value)
    return value.xreplace({f: sp.Rational(f) for f in value.atoms(sp.Float)})


def equal_poly(actual, expected):
    assert tuple(map(str, actual.gens)) == tuple(map(str, expected.gens))
    assert str(actual.domain) == str(expected.domain)
    assert sp.cancel(exact(actual.as_expr()) - exact(expected.as_expr())) == 0


CASES = [
    ('RR.frac_field(y)', '1.5*(x+z)*(x+1)/(y+1)', '3*(x+z)*(z+1)/(y+2)', ('x', 'z')),
    ('RR.frac_field(y)', '1.5*(x+z)/(y+1)', '3*(x+z)/(y+2)', ('x', 'z')),
    ('RR.frac_field(y)', '-1.5*(x+z)/(y+1)', '3*(x+z)/(y+2)', ('x', 'z')),
    ('RR.frac_field(y)', '1.5*(x+z)/(y+1)', '-3*(x+z)/(y+2)', ('x', 'z')),
    ('RR.frac_field(y)', '-1.5*(x+z)/(y+1)', '-3*(x+z)/(y+2)', ('x', 'z')),
    ('RR.frac_field(y)', '1.5*x/(y+1)', '3*x/(y+2)', ('x', 'z')),
    ('RR.frac_field(y)', '1.5*z/(y+1)', '3*z/(y+2)', ('x', 'z')),
    ('RR.frac_field(y)', 'x+y*z+1', 'x-y*z+2', ('x', 'z')),
    ('RR.frac_field(y)', '(x+z)/(1-y)', '(x+z)/(y-2)', ('x', 'z')),
    ('RR.frac_field(y)', '1.5', '3', ('x', 'z')),
    ('RR.frac_field(y,t)', '(x+z+w)*(x+y)/(y+t)', '(x+z+w)*(z+t)/(y+1)', ('x', 'z', 'w')),
    ('RR.frac_field(y,t)', '(x+y*z+t*w)*(x+1)/(y*t+1)', '(x+y*z+t*w)*(z+1)/(y+t)', ('x', 'z', 'w')),
    ('RR[y]', '1.5*(x+z)', '3*(x+z)', ('x', 'z')),
    ('RR[y]', '1.5*x', '3*x', ('x', 'z')),
    ('RR[y]', '-1.5*y*(x+z)', '3*y*(x+z)', ('x', 'z')),
    ('RR[y,t]', '(x+y*z+t*w)*(x+1)', '(x+y*z+t*w)*(z+1)', ('x', 'z', 'w')),
]


@pytest.mark.parametrize('domain,first,second,names', CASES)
def test_real_gcd_cofactors_lcm_and_cancel_match_reference(domain, first, second, names):
    answers = []
    for module in (sx, sp):
        K = evaluate(module, domain)
        gens = tuple(module.Symbol(n) for n in names)
        F, G = (module.Poly(evaluate(module, e), *gens, domain=K) for e in (first, second))
        h, p, q = F.cofactors(G)
        assert sp.cancel(exact(h.as_expr())*exact(p.as_expr()) - exact(F.as_expr())) == 0
        assert sp.cancel(exact(h.as_expr())*exact(q.as_expr()) - exact(G.as_expr())) == 0
        coefficient, a, b = F.cancel(G)
        assert sp.cancel(exact(coefficient)*exact(a.as_expr())*exact(G.as_expr())
                         - exact(b.as_expr())*exact(F.as_expr())) == 0
        included = F.cancel(G, include=True)
        answers.append((F.gcd(G), h, p, q, F.lcm(G), coefficient, a, b, *included))
    for index, (actual, expected) in enumerate(zip(*answers)):
        if index == 5:
            assert sp.cancel(exact(actual) - exact(expected)) == 0
        else:
            equal_poly(actual, expected)


@pytest.mark.parametrize('domain,first,second,names', CASES)
def test_expression_and_tuple_cancel_preserve_sparse_conventions(domain, first, second, names):
    answers = []
    for module in (sx, sp):
        K = evaluate(module, domain)
        gens = tuple(module.Symbol(n) for n in names)
        f, g = (evaluate(module, e) for e in (first, second))
        value = module.cancel(f/g, *gens, domain=K)
        triple = module.cancel((f, g), *gens, domain=K)
        assert len(triple) == 3
        coefficient, numerator, denominator = triple
        assert sp.cancel(exact(coefficient)*exact(numerator)*exact(g) - exact(denominator)*exact(f)) == 0
        answers.append((value, *triple))
    for actual, expected in zip(*answers):
        assert sp.cancel(exact(actual) - exact(expected)) == 0


@pytest.mark.parametrize('domain', ['RR.frac_field(y)', 'RR[y]'])
@pytest.mark.parametrize('first,second', [
    ('0', '3*(x+z)'), ('0', '-3*(x+z)'),
    ('3*(x+z)', '0'), ('-3*(x+z)', '0'), ('0', '0'),
])
def test_zero_associates_use_real_signs_and_keep_zero_gcd(domain, first, second):
    answers = []
    for module in (sx, sp):
        x, z = module.symbols('x z')
        K = evaluate(module, domain)
        f, g = (evaluate(module, e) for e in (first, second))
        F, G = (module.Poly(e, x, z, domain=K) for e in (f, g))
        if first == second == '0' and module is sp:
            with pytest.raises(ZeroDivisionError):
                F.lcm(G)
            lcm = module.Poly(0, x, z, domain=K)
        else:
            lcm = F.lcm(G)
        answers.append((F.gcd(G), *F.cofactors(G), lcm, F.cancel(G),
                        F.cancel(G, include=True), module.cancel((f, g), x, z, domain=K)))
    actual, expected = answers
    for a, b in zip(actual[:5], expected[:5]):
        equal_poly(a, b)
    assert sp.cancel(exact(actual[5][0]) - exact(expected[5][0])) == 0
    for a, b in zip(actual[5][1:] + actual[6], expected[5][1:] + expected[6]):
        equal_poly(a, b)
    assert tuple(map(exact, actual[7])) == tuple(map(exact, expected[7]))


@pytest.mark.parametrize('precision', [140, 211])
@pytest.mark.parametrize('ring', [False, True])
def test_high_precision_real_cofactors_reconstruct_actual_binary_coefficients(precision, ring):
    x, z, y = sx.symbols('x z y')
    base = sx.RealField(precision)
    K = base.poly_ring(y) if ring else base.frac_field(y)
    a = sx.Float('1.2345678901234567890123456789', precision=precision)
    c = a*y if ring else a/(y+1)
    F = sx.Poly(c*(x+z)*(x+1), x, z, domain=K)
    G = sx.Poly(2*c*(x+z)*(z+1), x, z, domain=K)
    h, p, q = F.cofactors(G)
    X, Z = sp.symbols('x z')
    for original, cofactor in ((F, p), (G, q)):
        source = sp.Poly(exact(original.as_expr()), X, Z)
        product = sp.Poly(exact(h.as_expr())*exact(cofactor.as_expr()), X, Z)
        assert set(product.monoms()) == set(source.monoms())
        for monomial, coefficient in source.terms():
            residual = sp.cancel((product.coeff_monomial(monomial) - coefficient)/coefficient)
            assert not residual.free_symbols
            assert abs(residual) < sp.Rational(1, 2)**(precision-4)
    for polynomial in (h, p, q, *F.cancel(G, include=True)):
        assert polynomial.domain.dom.precision == precision
        assert {f._prec for c in polynomial.coeffs() for f in c.atoms(sx.Float)} <= {dps_to_prec(prec_to_dps(precision))}
    result = sx.cancel(F.as_expr()/G.as_expr(), x, z, domain=K)
    assert sp.cancel(exact(result) - (X+1)/(2*(Z+1))) == 0


def test_parameter_identity_and_assumptions_survive_field_conversion():
    x, z = sx.symbols('x z')
    y = sx.Symbol('parameter', positive=True)
    K = sx.RealField(140).frac_field(y)
    F = sx.Poly((x+z)*(x+1)/(y+1), x, z, domain=K)
    G = sx.Poly((x+z)*(z+1)/(y+2), x, z, domain=K)
    for polynomial in (*F.cofactors(G), F.lcm(G), *F.cancel(G)[1:]):
        assert polynomial.domain.symbols == (y,)
        assert polynomial.gens == (x, z)
        assert polynomial.domain.dom.precision == 140


@pytest.mark.parametrize('source_domain', [None, 'QQ.frac_field(y)', 'QQ_I.frac_field(y)',
                                          'ZZ.frac_field(y)', 'QQ_I[y]', 'ZZ_I[y]'])
def test_cancel_polynomial_outputs_keep_the_reference_domains(source_domain):
    answers = []
    for module in (sx, sp):
        x, z, y = module.symbols('x z y')
        f, g = (x+z)*(x+1)/(y+1), (x+z)*(z+1)/(y+2)
        K = module.RR.frac_field(y)
        if source_domain:
            if source_domain.endswith(']'):
                f, g = (x+z)*(x+y), (x+z)*(z+y)
            original = tuple(module.Poly(e, x, z, domain=evaluate(module, source_domain)) for e in (f, g))
            result = module.cancel(original, domain=K)
        else:
            result = module.cancel((f, g), gens=(x, z), domain=K, polys=True)
        assert result[0] == 1
        assert all(isinstance(p, module.Poly) for p in result[1:])
        answers.append(result)
    for actual, expected in zip(answers[0][1:], answers[1][1:]):
        equal_poly(actual, expected)


@pytest.mark.parametrize('options', [{'polys': True}, {'field': True}, {'greedy': False}])
def test_cancel_field_option_errors_remain_explicit(options):
    for module in (sx, sp):
        x, z, y = module.symbols('x z y')
        with pytest.raises(module.polys.polyerrors.OptionError):
            module.cancel((x+z, x-z), x, z, domain=module.RR.frac_field(y), **options)


@pytest.mark.parametrize('precision', [53, 140])
@pytest.mark.parametrize('source', ['decimal', 'tiny', 'near_integer'])
@pytest.mark.parametrize('domain', ['QQ.frac_field(y)', 'QQ[y]', 'ZZ.frac_field(y)',
                                     'QQ_I.frac_field(y)', 'QQ_I[y]', 'ZZ_I.frac_field(y)'])
def test_explicit_exact_function_domain_conversion_uses_full_source_float(precision, source, domain):
    answers = []
    sources = []
    for module in (sx, sp):
        x, z, y = module.symbols('x z y')
        argument = ('0.1' if source == 'decimal' else '1e-60' if source == 'tiny'
                    else 2**100 + module.Rational(1, 2**20))
        value = module.Float(argument, precision=precision)
        sources.append(exact(value))
        K = evaluate(module, domain)
        result = module.Poly(value*x + y, x, z, domain=K)
        answers.append(result)
    equal_poly(*answers)
    assert sources[0] == sources[1]
    if source == 'decimal':
        assert exact(answers[0].LC()) == sp.Rational(1, 10)
    else:
        # Tiny values must not become zero. A 140-bit value 2**100+2**-20
        # must not be rounded to 53 bits before the denominator is bounded.
        assert exact(answers[0].LC()) == sources[0]
        if source == 'near_integer' and precision == 140:
            assert exact(answers[0].LC()) != 2**100


@pytest.mark.parametrize('source', ['Float(3)', 'Float("0.1")', 'Float("1e-60",precision=140)',
                                     'Float(1+Rational(1,2**100),precision=140)'])
def test_explicit_integer_ring_float_conversion_checks_the_unrounded_input(source):
    for module in (sx, sp):
        x, z, y = module.symbols('x z y')
        value = evaluate(module, source)
        if source == 'Float(3)':
            assert module.Poly(value*x, x, z, domain=module.ZZ[y]).as_expr() == 3*x
        else:
            with pytest.raises(module.polys.polyerrors.CoercionFailed):
                module.Poly(value*x, x, z, domain=module.ZZ[y])


@pytest.mark.parametrize('source', ['I + Float("0.1")', 'oo', 'nan'])
def test_exact_function_domains_still_reject_nonreal_or_nonfinite_values(source):
    for module in (sx, sp):
        x, z, y = module.symbols('x z y')
        value = evaluate(module, source)
        with pytest.raises(module.polys.polyerrors.CoercionFailed):
            module.Poly(value*x, x, z, domain=module.QQ.frac_field(y))


def test_float_exponents_are_not_rewritten_to_integer_powers_by_coefficient_conversion():
    for module in (sx, sp):
        x, z, y = module.symbols('x z y')
        with pytest.raises(module.polys.polyerrors.CoercionFailed):
            module.Poly(x*y**module.Float(2), x, z, domain=module.QQ.frac_field(y))


def test_default_exact_conversion_bound_does_not_enter_selected_140_bit_arithmetic():
    x, z, y = sx.symbols('x z y')
    value = sx.Float(1+sx.Rational(1, 2**100), precision=140)
    exact_poly = sx.Poly(value*x, x, z, domain=sx.QQ.frac_field(y))
    # This is the deliberate default-RR denominator bound of QQ.from_sympy.
    assert exact(exact_poly.LC()) == 1
    K = sx.RealField(140).frac_field(y)
    F = sx.Poly(value*(x+z), x, z, domain=K)
    G = sx.Poly(x+z, x, z, domain=K)
    _, p, q = F.cofactors(G)
    assert exact(p.LC()) == exact(value) == 1+sp.Rational(1, 2**100)
    assert exact(q.LC()) == 1
