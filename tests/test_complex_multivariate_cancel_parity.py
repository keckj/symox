"""Complete Gaussian exact-ground routes for multivariate inexact polynomials."""
import pytest
import sympy as sp
import symox as sx


def evaluate(module, source):
    namespace = vars(module).copy()
    namespace.update(zip(('x', 'z', 'w', 'y', 't'), module.symbols('x z w y t')))
    return eval(source, namespace)


def exact(expression):
    """The actual stored binary coefficients, without further rounding."""
    if isinstance(expression, sx.Basic):
        expression = sx.to_sympy(expression)
    expression = sp.sympify(expression)
    return expression.xreplace({f: sp.Rational(f) for f in expression.atoms(sp.Float)})


def assert_same_poly(actual, expected):
    assert tuple(map(str, actual.gens)) == tuple(map(str, expected.gens))
    assert str(actual.domain) == str(expected.domain)
    assert sp.cancel(exact(actual.as_expr()) - exact(expected.as_expr())) == 0


CASES = [
    ('CC', '(1 + I)*(x + z)*(x + 1)', '(2 - I)*(x + z)*(z + 2)', ('x', 'z')),
    ('CC', '1.5*(x + z)', '3*(x + z)', ('x', 'z')),
    ('CC', 'x + I*z + 1', 'x - I*z + 2', ('x', 'z')),
    ('CC', '(x + z + w)*(x + 1)', '(x + z + w)*(z + 1)', ('x', 'z', 'w')),
    ('CC[y]', '1.5*(x + z)', '3*(x + z)', ('x', 'z')),
    ('CC[y]', '1.5*x', '3*x', ('x', 'z')),
    ('CC[y]', '1.5*z', '3*z', ('x', 'z')),
    ('CC[y]', '1.5*y', '3*y**2', ('x', 'z')),
    ('CC[y]', '1.5*(z + 1)', '3*(z + 1)', ('x', 'z')),
    ('CC[y]', '(1 + I)*(x + z)', '2*(x + z)', ('x', 'z')),
    ('CC[y]', '1.5*(x + y*z + 1)*(x + z)', '3*(x + y*z + 1)*(z + y)', ('x', 'z')),
    ('CC[y, t]', '(x + y*z + t*w)*(x + I*z)', '(x + y*z + t*w)*(z + w)', ('x', 'z', 'w')),
    ('CC.frac_field(y)', '1.5*(x + z)', '3*(x + z)', ('x', 'z')),
    ('CC.frac_field(y)', '(x + z)*(x + 1)/(y + 1)', '(x + z)*(z + 1)/(y + 2)', ('x', 'z')),
    ('CC.frac_field(y, t)', '(x + z + w)*(x + y)/(y + t)', '(x + z + w)*(z + t)/(y + 1)', ('x', 'z', 'w')),
]


@pytest.mark.parametrize('domain, first, second, names', CASES)
def test_gcd_cofactors_cancel_and_lcm_match_the_reference(domain, first, second, names):
    results = []
    for module in (sx, sp):
        K = evaluate(module, domain)
        gens = tuple(module.Symbol(n) for n in names)
        F = module.Poly(evaluate(module, first), *gens, domain=K)
        G = module.Poly(evaluate(module, second), *gens, domain=K)
        h, p, q = F.cofactors(G)
        # Prove both identities over exact rationals, not a numeric sample.
        assert sp.cancel(exact(h.as_expr())*exact(p.as_expr()) - exact(F.as_expr())) == 0
        assert sp.cancel(exact(h.as_expr())*exact(q.as_expr()) - exact(G.as_expr())) == 0
        cancelled = F.cancel(G)
        c, a, b = cancelled
        assert sp.cancel(exact(c)*exact(a.as_expr())*exact(G.as_expr())
                         - exact(b.as_expr())*exact(F.as_expr())) == 0
        included = F.cancel(G, include=True)
        assert sp.cancel(exact(included[0].as_expr())*exact(G.as_expr())
                         - exact(included[1].as_expr())*exact(F.as_expr())) == 0
        results.append((F.gcd(G), (h, p, q), F.lcm(G), cancelled, included))
    actual, expected = results
    assert_same_poly(actual[0], expected[0])
    for a, b in zip(actual[1], expected[1]):
        assert_same_poly(a, b)
    assert_same_poly(actual[2], expected[2])
    assert sp.cancel(exact(actual[3][0]) - exact(expected[3][0])) == 0
    for a, b in zip(actual[3][1:], expected[3][1:]):
        assert_same_poly(a, b)
    for a, b in zip(actual[4], expected[4]):
        assert_same_poly(a, b)


@pytest.mark.parametrize('domain', ['CC', 'CC[y]', 'CC.frac_field(y)'])
@pytest.mark.parametrize('first,second', [
    ('0', '(1 + I)*(x + z)'),
    ('(1 + I)*(x + z)', '0'),
    ('0', '0'),
    ('1.5', '3'),
])
def test_zero_and_constant_polynomials_do_not_fake_a_unit_gcd(domain, first, second):
    polynomials = []
    for module in (sx, sp):
        x, z = module.symbols('x z')
        K = evaluate(module, domain)
        F, G = (module.Poly(evaluate(module, e), x, z, domain=K) for e in (first, second))
        if module is sp and first == second == '0':
            # Reference dmp_rr_lcm/dmp_ff_lcm divides 0 by gcd(0, 0).
            # Symox deliberately returns the mathematical zero LCM.
            with pytest.raises(ZeroDivisionError):
                F.lcm(G)
            least_common_multiple = module.Poly(0, x, z, domain=K)
        else:
            least_common_multiple = F.lcm(G)
        polynomials.append((F.gcd(G), least_common_multiple, F.cofactors(G), F.cancel(G)))
    actual, expected = polynomials
    assert_same_poly(actual[0], expected[0])
    assert_same_poly(actual[1], expected[1])
    for a, b in zip(actual[2], expected[2]):
        assert_same_poly(a, b)
    assert exact(actual[3][0]) == exact(expected[3][0])
    for a, b in zip(actual[3][1:], expected[3][1:]):
        assert_same_poly(a, b)


@pytest.mark.parametrize('domain,first,second,names', CASES)
def test_public_expression_cancel_and_tuple_protocol(domain, first, second, names):
    answers = []
    for module in (sx, sp):
        K = evaluate(module, domain)
        gens = tuple(module.Symbol(n) for n in names)
        f, g = evaluate(module, first), evaluate(module, second)
        value = module.cancel(f/g, *gens, domain=K)
        coefficient, numerator, denominator = module.cancel((f, g), *gens, domain=K)
        assert sp.cancel(exact(coefficient)*exact(numerator)*exact(g) - exact(denominator)*exact(f)) == 0
        answers.append((value, coefficient, numerator, denominator))
    for actual, expected in zip(*answers):
        assert sp.cancel(exact(actual) - exact(expected)) == 0


@pytest.mark.parametrize('convert', [False, True])
def test_multivariate_fraction_field_clears_parameter_denominators(convert):
    answers = []
    for module in (sx, sp):
        x, z, y = module.symbols('x z y')
        K = module.CC.frac_field(y)
        F = module.Poly((x + z)/(y + 1) + (x*z + 1)/(y + 2), x, z, domain=K)
        coefficient, result = F.clear_denoms(convert=convert)
        assert sp.cancel(exact(coefficient)*exact(F.as_expr()) - exact(result.as_expr())) == 0
        answers.append((coefficient, result))
    assert sp.cancel(exact(answers[0][0]) - exact(answers[1][0])) == 0
    assert_same_poly(answers[0][1], answers[1][1])


@pytest.mark.parametrize('kind', ['scalar', 'ring', 'field'])
def test_140_bit_complex_coefficients_are_preserved_in_multivariate_cofactors(kind):
    x, z, y = sx.symbols('x z y')
    K = sx.ComplexField(140)
    if kind == 'ring':
        K = K.poly_ring(y)
    elif kind == 'field':
        K = K.frac_field(y)
    a = sx.Float('1.2345678901234567890123456789', 40)*(1 + sx.I)
    coefficient = a*y if kind == 'ring' else a/(y + 1) if kind == 'field' else a
    F = sx.Poly(coefficient*(x + z)*(x + 1), x, z, domain=K)
    G = sx.Poly(2*coefficient*(x + z)*(z + 1), x, z, domain=K)
    h, p, q = F.cofactors(G)
    X, Z = sp.symbols('x z')
    for original, cofactor in ((F, p), (G, q)):
        source = sp.Poly(exact(original.as_expr()), X, Z)
        product = sp.Poly(exact(h.as_expr())*exact(cofactor.as_expr()), X, Z)
        assert set(product.monoms()) == set(source.monoms())
        for monomial, stored in source.terms():
            # Separate 140-bit rounding of h and its cofactor need not
            # reconstruct an identical rational. Bound each relative error
            # using exact binary inputs, allowing 16 ulps for normalization.
            residual = sp.cancel((product.coeff_monomial(monomial) - stored)/stored)
            assert not residual.free_symbols
            squared_error = sp.cancel(sp.re(residual)**2 + sp.im(residual)**2)
            assert squared_error < sp.Rational(1, 2)**272
    for polynomial in (h, p, q, *F.cancel(G, include=True)):
        ground = polynomial.domain if kind == 'scalar' else polynomial.domain.dom
        assert ground.precision == 140
        assert all(f._prec == 140 for c in polynomial.coeffs() for f in c.atoms(sx.Float))
    result = sx.cancel(F.as_expr()/G.as_expr(), x, z, domain=K)
    expected = (sp.Symbol('x') + 1)/(2*(sp.Symbol('z') + 1))
    assert sp.cancel(exact(result) - expected) == 0


def test_parameter_denominator_outside_ring_keeps_coercion_failure():
    for module in (sx, sp):
        x, z, y = module.symbols('x z y')
        with pytest.raises(module.polys.polyerrors.CoercionFailed):
            module.Poly((x + z)/(y + 1), x, z, domain=module.CC[y])


@pytest.mark.parametrize('domain', ['CC', 'CC[y]', 'CC.frac_field(y)'])
@pytest.mark.parametrize('first,second', [('0', '2*(x+z)'), ('0', '0'), ('x+z', '0')])
def test_expression_cancel_zero_protocol(domain, first, second):
    answers = []
    for module in (sx, sp):
        x, z = module.symbols('x z')
        answers.append(module.cancel((evaluate(module, first), evaluate(module, second)),
                                    x, z, domain=evaluate(module, domain)))
    assert tuple(map(exact, answers[0])) == tuple(map(exact, answers[1]))


@pytest.mark.parametrize('ground', ['CC', 'RR', 'ComplexField(140)', 'RealField(140)'])
def test_inexact_composite_domains_have_their_associated_ring(ground):
    for module in (sx, sp):
        y = module.Symbol('y')
        base = evaluate(module, ground)
        assert not base.has_assoc_Ring
        for domain in (base.poly_ring(y), base.frac_field(y)):
            assert domain.has_assoc_Ring
            assert domain.get_ring() == base.poly_ring(y)


@pytest.mark.parametrize('domain', ['ZZ', 'QQ', 'ZZ_I', 'QQ_I', 'RR', 'RR[y]', 'RR.frac_field(y)'])
@pytest.mark.parametrize('multivariate', [False, True])
def test_sparse_adapter_preserves_other_domain_cofactor_routes(domain, multivariate):
    x, z = sx.symbols('x z')
    gens = (x, z) if multivariate else (x,)
    shared = x + z if multivariate else x + 1
    K = evaluate(sx, domain)
    F, G = (sx.Poly(e, *gens, domain=K) for e in (shared*(x + 2), 2*shared))
    assert F._sparse_cofactors(G) == F.cofactors(G)
    with pytest.raises(TypeError):
        F._sparse_cofactors(None)
