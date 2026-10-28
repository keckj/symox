"""Large total degree is legal when each stored coordinate remains legal."""
import pytest
import sympy as sp
import symox as sx

N = 1 << 31


def exact_dict(p):
    result = {}
    for powers, coefficient in p.as_dict().items():
        expression = sx.to_sympy(coefficient)
        result[powers] = expression.xreplace({f: sp.Rational(f) for f in expression.atoms(sp.Float)})
    return result


def domain(module, name, y):
    return {
        'ZZ': lambda: module.ZZ,
        'QQ': lambda: module.QQ,
        'QQ_I': lambda: module.QQ_I,
        'QQ(y)': lambda: module.QQ.frac_field(y),
        'RR140': lambda: module.RealField(140),
        'CC140': lambda: module.ComplexField(140),
        'EX': lambda: module.EX,
    }[name]()


@pytest.mark.parametrize('name', ['ZZ', 'QQ', 'QQ_I', 'QQ(y)', 'RR140'])
@pytest.mark.parametrize('support', [{(1, 1): 1}, {(1, 1): 1, (1, 0): 2}, {(1, 1): 1, (0, 0): 3}])
def test_total_degree_and_homogeneity_follow_exact_scaled_support(name, support):
    x, y, a = sx.symbols('x y a')
    sparse = {tuple(N*k for k in powers): value for powers, value in support.items()}
    p = sx.Poly.from_dict(sparse, x, y, domain=domain(sx, name, a))
    # Scaling the reference's small sparse support proves the large degree
    # without asking its dense DMP representation to allocate N zero slots.
    X, Y = sp.symbols('x y')
    reference = sp.Poly.from_dict(support, X, Y)
    assert p.total_degree() == N*reference.total_degree()
    assert sx.total_degree(p) == p.total_degree()
    assert p.degree_list() == (max(e[0] for e in sparse), max(e[1] for e in sparse))
    assert p.is_homogeneous == reference.is_homogeneous
    assert not p.is_linear and not p.is_quadratic
    expected_order = N*reference.homogeneous_order() if reference.is_homogeneous else None
    assert p.homogeneous_order() == expected_order
    assert exact_dict(p) == sparse


@pytest.mark.parametrize('name', ['ZZ', 'QQ', 'QQ_I', 'QQ(y)', 'RR140'])
def test_homogenization_preserves_representable_coordinates_and_coefficients(name):
    x, y, z, a = sx.symbols('x y z a')
    K = domain(sx, name, a)
    p = sx.Poly.from_dict({(N, N): 1, (N, 0): 2}, x, y, domain=K)
    result = p.homogenize(z)
    assert exact_dict(result) == {(N, N, 0): 1, (N, 0, N): 2}
    assert result.gens == (x, y, z) and result.domain == p.domain
    assert result.total_degree() == 2*N
    assert result.is_homogeneous and result.homogeneous_order() == 2*N
    assert result.eval(z, 1) == p
    # Homogenizing in an existing generator can combine terms exactly.
    assert exact_dict(p.homogenize(y)) == {(N, N): 3}
    assert p.homogenize(y).domain == p.domain
    with pytest.raises(NotImplementedError, match='exponent storage limit'):
        p.homogenize(x)
    q = sx.Poly.from_dict({(N, N): 1, (0, 0): 1}, x, y, domain=K)
    with pytest.raises(NotImplementedError, match='exponent storage limit'):
        q.homogenize(z)


@pytest.mark.parametrize('order', ['grlex', 'grevlex'])
def test_domain_element_display_uses_the_wide_total_degree(order):
    x, y = sx.symbols('x y')
    X, Y = sp.symbols('x y')
    K = sx.QQ.poly_ring(x, y, order=order)
    reference = sp.QQ.poly_ring(X, Y, order=order)
    expression = x**(2*N-1) + x**N*y**N
    expected = X**(2*N-1) + X**N*Y**N
    assert str(K.from_sympy(expression)) == str(reference.from_sympy(expected))
    assert sx.to_sympy(K.from_sympy(expression).as_expr()) == expected


def test_three_variable_wide_degree_and_zero_conventions():
    x, y, z = sx.symbols('x y z')
    p = sx.Poly.from_dict({(2*N-1, 2*N-1, 2*N-1): 1}, x, y, z)
    assert p.total_degree() == 3*(2*N-1)
    assert p.homogeneous_order() == p.total_degree()
    assert p.homogenize(z) == p
    assert sx.Poly(0, x, y).total_degree() == 0
    assert sx.Poly(0, x, y).homogeneous_order() == -sx.oo


@pytest.mark.parametrize('name', ['CC140', 'EX'])
def test_dense_expression_domains_keep_regular_degree_and_homogenization(name):
    # These grounds use a dense DMP representation; giant coordinates are
    # outside this test's allocation budget, independently of total degrees.
    x, y, z, a = sx.symbols('x y z a')
    p = sx.Poly.from_dict({(3, 3): 1, (3, 0): 2}, x, y, domain=domain(sx, name, a))
    assert p.total_degree() == 6 and sx.total_degree(p) == 6
    result = p.homogenize(z)
    assert exact_dict(result) == {(3, 3, 0): 1, (3, 0, 3): 2}
    assert result.domain == p.domain
    assert result.homogeneous_order() == 6
    if name == 'CC140':
        assert {f._prec for c in result.coeffs() for f in c.atoms(sx.Float)} == {140}


@pytest.mark.parametrize('complex_', [False, True])
def test_homogenization_keeps_actual_140_bit_coefficient_values(complex_):
    x, y, z = sx.symbols('x y z')
    K = sx.ComplexField(140) if complex_ else sx.RealField(140)
    value = sx.Float('1.23456789012345678901234567890123456789', precision=140)
    value *= sx.I if complex_ else 1
    p = sx.Poly.from_dict({(3, 3): value, (3, 0): 2}, x, y, domain=K)
    result = p.homogenize(z)
    expected = sx.to_sympy(value)
    expected = expected.xreplace({f: sp.Rational(f) for f in expected.atoms(sp.Float)})
    assert exact_dict(result) == {(3, 3, 0): expected, (3, 0, 3): 2}
    assert result.domain == p.domain
    assert {f._prec for c in result.coeffs() for f in c.atoms(sx.Float)} == {140}


@pytest.mark.parametrize('gens', [(), ('x',), ('y',), ('x', 'y'), ('y', 'x')])
def test_total_degree_keeps_the_expression_route_for_different_generators(gens):
    x, y = sx.symbols('x y')
    X, Y = sp.symbols('x y')
    p = sx.Poly(x**2*y + x*y**2, x, y)
    reference = sp.Poly(X**2*Y + X*Y**2, X, Y)
    selected = {'x': x, 'y': y}
    oracle = {'x': X, 'y': Y}
    assert sx.total_degree(p, *(selected[g] for g in gens)) == sp.total_degree(reference, *(oracle[g] for g in gens))
