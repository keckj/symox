"""Composite coefficient domains follow the reference inference and algorithms."""
import pytest
import symox as sm


def both(case):
    sy = pytest.importorskip('sympy')
    return case(sm), case(sy)


def same(case):
    ours, reference = both(case)
    assert str(ours) == str(reference)


@pytest.mark.parametrize('name', ['gcdex', 'half_gcdex', 'invert'])
def test_extended_euclid_over_rational_function_fields(name):
    def case(S):
        x, a = S.symbols('x a')
        f = S.Poly(x**2 + a, x, domain='QQ(a)')
        g = S.Poly(x + a, x, domain='QQ(a)')
        return g.invert(f) if name == 'invert' else getattr(f, name)(g)
    same(case)


def test_extended_euclid_errors_match_the_reference():
    sy = pytest.importorskip('sympy')
    for S in (sm, sy):
        x, a = S.symbols('x a')
        g = S.Poly(x + a, x, domain='QQ(a)')
        with pytest.raises(S.polys.polyerrors.NotInvertible):
            g.invert(S.Poly(x**2 - a**2, x, domain='QQ(a)'))
        f, g = S.Poly(x + a, x, domain='ZZ[a]'), S.Poly(x, x, domain='ZZ[a]')
        with pytest.raises(S.polys.polyerrors.DomainError, match='ground domain must be a field'):
            f.gcdex(g, auto=False)
        with pytest.raises(S.polys.polyerrors.DomainError, match=r'half extended GCD over ZZ\[a\]'):
            f.half_gcdex(g, auto=False)
        with pytest.raises(ValueError, match='univariate polynomial expected'):
            S.gcdex(x + a, x)


@pytest.mark.parametrize('call', [
    lambda S, x, a: S.gcdex(S.Poly(x + a, x, domain='ZZ(a)'), S.Poly(x, x, domain='ZZ(a)')),
    lambda S, x, a: S.gcdex(x + a, x, x),
    lambda S, x, a: S.gcdex(S.Poly(x + 1, x), x**2, x),
    lambda S, x, a: S.half_gcdex(x**2 + a, x + 1, x),
    lambda S, x, a: S.invert(x + a, x**2, x),
    lambda S, x, a: S.gcdex(S.pi, 2),
    lambda S, x, a: S.gcdex(2, 4),
    lambda S, x, a: S.half_gcdex(S.Rational(1, 2), 3),
    lambda S, x, a: S.ZZ.to_sympy(S.ZZ.invert(S.ZZ(1), S.ZZ(-3))),
    lambda S, x, a: [S.ZZ.to_sympy(c) for c in S.ZZ.gcdex(S.ZZ(0), S.ZZ(0))],
    lambda S, x, a: [S.RR.to_sympy(c) for c in S.RR.gcdex(S.RR(1), S.RR(2))],
    lambda S, x, a: S.QQ_I.invert(S.QQ_I(1), S.QQ_I.convert(1 + S.I)),
    lambda S, x, a: S.polys.constructor.construct_domain([S.pi / 2, 1])[0],
    lambda S, x, a: S.polys.constructor.construct_domain([S.pi, 2], composite=False)[0],
])
def test_functional_euclid_and_ground_domains(call):
    def case(S):
        x, a = S.symbols('x a')
        try:
            return call(S, x, a)
        except Exception as e:
            return type(e).__name__
    same(case)


@pytest.mark.parametrize('include', [True, False])
def test_cancel_over_fraction_fields(include):
    def case(S):
        x, a = S.symbols('x a')
        p = S.Poly(x**2 / a - a, x, domain='ZZ(a)')
        q = S.Poly(x / (2 * a) + 1, x, domain='ZZ(a)')
        return p.cancel(q, include=include)
    same(case)


@pytest.mark.parametrize('source', [
    lambda S, x, a: (S.Poly(x + S.pi, x), S.Poly(x / S.pi, x)),
    lambda S, x, a: (S.Poly(x + S.exp(-a), x), S.Poly(x * S.exp(a), x)),
    lambda S, x, a: (S.Poly(x * a, x), S.Poly(x / S.sqrt(a), x)),
    lambda S, x, a: (S.Poly(x / S.sqrt(a), x), S.Poly(x * a, x)),
])
def test_unified_generators_keep_receiver_order(source):
    def case(S):
        x, a = S.symbols('x a')
        f, g = source(S, x, a)
        return (f * g).domain, (g * f).domain
    same(case)


def test_groebner_infers_transcendental_fraction_fields():
    def case(S):
        x, y = S.symbols('x y')
        return S.groebner([x * y - 1 / S.pi, x - S.pi * y], x, y)
    same(case)


@pytest.mark.parametrize('source', [
    'a*x**4 - 16*a**5',
    'x**2 + (a + 1/a)*x + 1',
    'x**4 - a**2*x**2 - 2*a**2 - 1',
    'a*x**2 + b*x + c',
    '(x - a)**2*(x + b)',
])
def test_symbolic_roots_and_solve(source):
    def case(S):
        x = S.Symbol('x')
        e = S.sympify(source)
        return sorted(map(str, S.roots(e, x).items())), sorted(map(str, S.solve(e, x)))
    same(case)


@pytest.mark.parametrize('function', ['factor', 'sqf'])
@pytest.mark.parametrize('source', ['x**2 - 1/pi**2', 'a*x**2 - a', 'x**3 - 2*x**2 + x'])
def test_factoring_a_poly_returns_a_poly_with_inferred_domain(function, source):
    def case(S):
        x = S.Symbol('x')
        r = getattr(S, function)(S.Poly(S.sympify(source), x))
        return type(r).__name__, r
    same(case)
