"""Nested polynomial rings, extension=True resultants and domain conversions follow the reference."""
import pytest
import symox as sm


def namespace(S):
    g = {k: getattr(S, k) for k in dir(S) if not k.startswith('_')}
    d = S.polys.domains
    g.update({k: getattr(d, k) for k in ['ZZ', 'QQ', 'ZZ_I', 'QQ_I', 'RR', 'CC', 'EX', 'EXRAW']})
    exec("x, y, z, t = symbols('x y z t')", g)
    return g


def outcome(S, case):
    try:
        r = eval(case, namespace(S))
    except Exception as e:
        return 'EXC', type(e).__name__
    return str(r), repr(r)


NESTED_CASES = [
    'ZZ[x][y]',
    "PolynomialRing(ZZ[x], 'y')",
    'repr(ZZ[x][y])',
    'ZZ[x][y].domain',
    'ZZ[x][y].gens',
    'ZZ[x][y].symbols',
    'ZZ_I[x]',
    'QQ_I[x]',
    'QQ.algebraic_field(sqrt(2))[x]',
    'QQ.algebraic_field(sqrt(2)).poly_ring(x, y)',
    'ZZ_I.frac_field(x)',
    'QQ.algebraic_field(sqrt(2)).frac_field(x)',
    'ZZ[x].frac_field(y)',
    'ZZ[x][y](x*y+1)',
    'ZZ[x][y](x*y+x+1)',
    'repr(ZZ[x][y](x*y+x+1))',
    'ZZ[x][y](x*y+x+1)**2',
    'ZZ_I[x](x**2 + I)',
    'ZZ_I[x](x**2 + I)*ZZ_I[x](x - I)',
    'QQ_I[x](x/2 + I/3)',
    'QQ.algebraic_field(sqrt(2))[x](sqrt(2)*x + 1)',
    'QQ.algebraic_field(sqrt(2))[x](sqrt(2)*x + 1)**2',
    'ZZ[x][y].to_sympy(ZZ[x][y](x*y+x+1))',
    'ZZ[x][y].from_sympy(x*y+x+1)',
    'ZZ[x][y].convert(x*y)',
    'ZZ[x][y].unify(ZZ[x])',
    'ZZ[x][y].unify(QQ)',
    'ZZ[x][y].unify(ZZ[y])',
    'ZZ[x][y].unify(ZZ[x,y])',
    'ZZ_I[x].unify(QQ)',
    'ZZ_I[x].unify(QQ_I)',
    'ZZ_I[x].unify(ZZ[x])',
    'ZZ_I[x].unify(ZZ[y])',
    'QQ.algebraic_field(sqrt(2))[x].unify(QQ[x])',
    'QQ.algebraic_field(sqrt(2))[x].unify(ZZ_I)',
    'ZZ[x][y].get_field()',
    'ZZ_I[x].get_field()',
    'ZZ[x][y].get_ring()',
    'ZZ[x][y].is_PolynomialRing',
    'ZZ[x][y].ngens',
    'ZZ[x][y] == PolynomialRing(ZZ[x], y)',
    'ZZ[x][y] == ZZ[x,y]',
    'ZZ[x][y].dom',
    'ZZ_I[x].domain',
    'Poly(x*y+x+1, y, domain=ZZ[x])',
    'Poly(x*y*z+x+1, z, domain=ZZ[x][y])',
    'Poly(x**2 + I, x, domain=ZZ_I)',
    'Poly(x*y*z + 1, z, domain=ZZ[x][y]).rep',
    'Poly(x*y*z + 1, z, domain=ZZ[x][y]).as_expr()',
    'Poly(x*y*z + 1, z, domain=ZZ[x][y]).coeffs()',
    'Poly(x**2*y**2 - 1, y, domain=ZZ[x][z]).factor_list()',
    'ZZ[x][y].is_Field',
    'ZZ_I[x].characteristic()',
    'ZZ[x][y].is_Exact',
    'ZZ[x][y](x*y) // ZZ[x][y](y)',
    'ZZ[x][y].exquo(ZZ[x][y](x*y), ZZ[x][y](x))',
    'ZZ[x][y].gcd(ZZ[x][y](x*y), ZZ[x][y](x**2))',
    'ZZ_I[x].gcd(ZZ_I[x](x**2+1), ZZ_I[x](x+I))',
    'ZZ[x][y].is_unit(ZZ[x][y](-1))',
    'ZZ_I[x](x+I) == ZZ_I[x](x+I)',
    'ZZ_I[x].one',
    'ZZ_I[x].zero',
    'ZZ[x][y].gens[0]',
    'ZZ_I[x](I*x).as_expr()',
    'ZZ[x][y](2*x*y)',
    'ZZ[x][y](-x*y-1)',
    'ZZ[x][y](-x**2*y-x)',
    'QQ[x][y](x*y/2+y/3)',
    'ZZ[x][y].zero',
    'ZZ[x][y](y**2*(x+1)-x)',
    'ZZ[x][y][z](x*y*z+y+1)',
    'ZZ_I[x](-x)',
    'ZZ_I[x](2*I*x**2)',
    'Poly(x*y*z+x+1, z, domain=ZZ[x][y]).domain',
    'Poly(x*y*z + 1, z, domain=ZZ[x][y])**2',
    'Poly(x*y*z + 1, z, domain=ZZ[x][y]) + Poly(x*z, z, domain=ZZ[x][y])',
    'Poly(x*y*z + 1, z, domain=ZZ[x][y]) * 2',
    'Poly(x*y*z + 1, z, domain=ZZ[x][y]) + 1',
    'Poly(x**2*z**2 - 1, z, domain=ZZ[x][y]).factor_list()',
    'gcd(Poly(x*y*z**2 - x*y, z, domain=ZZ[x][y]), Poly(x*z - x, z, domain=ZZ[x][y]))',
    'div(Poly(x*z**2 + 1, z, domain=ZZ[x][y]), Poly(y*z + 1, z, domain=ZZ[x][y]))',
    'div(Poly(x*z**2 + 1, z, domain=QQ[x][y]), Poly(2*z + 1, z, domain=QQ[x][y]))',
    'Poly(x*z/2, z, domain=ZZ[x][y])',
    'Poly(x*z/2, z, domain=QQ[x][y])',
    'Poly(z/y, z, domain=ZZ[x][y])',
    'Poly(x*y*z + 1, z, domain=ZZ[x][y]).LC()',
    'Poly(x*y*z + 1, z, domain=ZZ[x][y]).eval(2)',
    'Poly(x*y*z + 1, z, domain=ZZ[x][y]).diff(z)',
    'Poly(x*y*z + 1, z, domain=ZZ[x][y]).get_domain()',
    'Poly(x*y*z + 1, z, domain=ZZ[x][y]).to_field()',
    'Poly(x*y*z + y, z, domain=ZZ[x][y]).content()',
    'Poly(x*y*z + y, z, domain=ZZ[x][y]).primitive()',
    'Poly(x*y*z + 1, z, domain=ZZ[x][y]).degree()',
    'resultant(Poly(x*z + y, z, domain=ZZ[x][y]), Poly(z**2 - x, z, domain=ZZ[x][y]))',
    'Poly(x*y*z+1, z, domain=ZZ[x][y]) == Poly(x*y*z+1, z, domain=ZZ[x,y])',
    'Poly(x*y*z+1, z, domain=ZZ[x][y]).unify(Poly(z, z, domain=QQ))',
    'sqf_list(Poly(x**2*z**2 + 2*x*z + 1, z, domain=ZZ[x][y]))',
    'Poly(x*y*z + 1, z, domain=ZZ[x][y]).to_field().domain',
    'Poly(x*y*z + 1, z, domain=ZZ[x][y]) * Poly(z, z, domain=ZZ)',
    'Poly(x*y*z + 1, z, domain=ZZ[x][y]) - Poly(z/2, z, domain=QQ)',
    'Poly(x*z + 1, z, domain=ZZ[x].frac_field(y))',
    'Poly(x*z/y + 1, z, domain=ZZ[x].frac_field(y))',
    'Poly(x*z/y + 1, z, domain=ZZ[x].frac_field(y)).domain',
]


@pytest.mark.parametrize('case', NESTED_CASES)
def test_nested_ring_case_matches_reference(case):
    sy = pytest.importorskip('sympy')
    assert outcome(sm, case) == outcome(sy, case)


def test_sage_style_polynomial_ring_conversion(monkeypatch):
    sy = pytest.importorskip('sympy')
    for S in (sm, sy):
        cls = __import__(S.__name__ + '.polys.domains.polynomialring', fromlist=['PolynomialRing']).PolynomialRing
        monkeypatch.setattr(cls, '_sage_', lambda self: (str(self.domain), ','.join(map(str, self.gens))), raising=False)
        x, y = S.symbols('x y')
        ZZx = S.ZZ[x]
        R = S.PolynomialRing(ZZx, 'y')
        assert str(R) == 'ZZ[x][y]'
        assert R.is_PolynomialRing
        assert str(R.domain) == 'ZZ[x]'
        assert ','.join(map(str, R.gens)) == 'y'
        assert ','.join(map(str, R.domain.gens)) == 'x'
        assert isinstance(R, cls) and isinstance(R.domain, cls)
        assert S.PolynomialRing(ZZx, 'y')._sage_() == ('ZZ[x]', 'y')


EXTENSION_CASES = [
    'resultant(x**2 - sqrt(2)*y, x - y, x, extension=True)',
    'resultant(x**2 - I*y, x - y, x, extension=True)',
    'resultant(x**2 - sqrt(2), x - y, x, extension=True)',
    'resultant(x**2 - sqrt(2)*y, x - y, x, extension=True, polys=True)',
    'subresultants(x**2 - sqrt(2)*y, x - y, x, extension=True)',
    'discriminant(x**2 - sqrt(2)*y*x + 1, x, extension=True)',
    'discriminant(x**3 - sqrt(2)*y, x, extension=True)',
    'gcd(x**2 - 2*y**2, x - sqrt(2)*y, x, extension=True)',
    'lcm(x**2 - 2*y**2, x - sqrt(2)*y, x, extension=True)',
    'cancel((x**2 - 2*y**2)/(x - sqrt(2)*y), extension=True)',
    '[p.domain for p in parallel_poly_from_expr([x**2 - sqrt(2)*y, x - y], x, extension=True)[0]]',
    '[p.domain for p in parallel_poly_from_expr([x**2 - sqrt(2), x - 1], x, extension=True)[0]]',
    'resultant(Poly(x**2 - sqrt(2)*y, x, extension=True), Poly(x - y, x))',
]


@pytest.mark.parametrize('case', EXTENSION_CASES)
def test_extension_true_with_parameters_matches_reference(case):
    sy = pytest.importorskip('sympy')
    assert outcome(sm, case) == outcome(sy, case)


@pytest.mark.parametrize('a, b', [(4, 0), (3, 0), (4, 7), (0.5, 3)])
def test_real_field_invert_matches_reference(a, b):
    sy = pytest.importorskip('sympy')
    ours = sm.polys.domains.RR.invert(a, b)
    reference = sy.polys.domains.RR.invert(a, b)
    assert float(ours) == float(reference)


def test_real_field_invert_errors_match_reference():
    sy = pytest.importorskip('sympy')
    for S in (sm, sy):
        RR = S.polys.domains.RR
        with pytest.raises(ZeroDivisionError):
            RR.invert(-2, 0)
        with pytest.raises(S.polys.polyerrors.NotInvertible):
            RR.invert(0, 0)


@pytest.mark.parametrize('domain', ['ZZ', 'QQ', 'RR', 'CC', 'EX', 'EXRAW', 'ZZ_I', 'QQ_I'])
@pytest.mark.parametrize('value', ['x+1', None, b'x', [1], (1,)])
def test_non_coefficients_are_rejected_like_reference(domain, value):
    sy = pytest.importorskip('sympy')
    for S in (sm, sy):
        with pytest.raises(S.polys.polyerrors.CoercionFailed):
            getattr(S.polys.domains, domain).convert(value)


def test_raw_expression_domain_conversions_match_reference():
    sy = pytest.importorskip('sympy')
    for S in (sm, sy):
        EXRAW = S.polys.domains.EXRAW
        assert EXRAW.to_sympy(EXRAW.convert(S.oo)) == S.oo
        assert EXRAW.to_sympy(EXRAW.convert(S.Symbol('x') + 1)) == S.Symbol('x') + 1
        with pytest.raises(S.polys.polyerrors.CoercionFailed, match='Expecting an Expr instance'):
            EXRAW.from_sympy('x')
        with pytest.raises(S.polys.polyerrors.CoercionFailed, match='is not in any domain'):
            S.polys.domains.QQ.convert(S.oo)
