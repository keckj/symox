"""Scalar coefficient calculus with noncommutative matrix factors."""
import pytest
import sympy as sp
import symox as sx


def expressions():
    x = sx.Symbol('x')
    a, b = sx.MatrixSymbol('A', 2, 2), sx.MatrixSymbol('B', 2, 2)
    return x, a, b, [a, a*b, 2*x*a+b, sx.exp(x)*a, (x*x+1)*a*b,
                     sx.sin(x)*a+sx.cos(x)*b, sx.Rational(2, 3)*a]


@pytest.mark.parametrize('index', range(7))
@pytest.mark.parametrize('method', ['as_terms', 'as_ordered_terms', 'as_coefficients_dict'])
def test_term_views(index, method):
    from test_matrixexpr_protocols import convert
    _, _, _, es = expressions()
    e = es[index]
    actual, expected = getattr(e, method)(), getattr(e.to_sympy(), method)()
    if method == 'as_terms':
        terms, generators = actual
        ref_terms, ref_generators = expected
        assert convert(generators) == ref_generators
        assert dict(convert(terms)) == dict(ref_terms)
    else:
        assert convert(actual) == expected


@pytest.mark.parametrize('index', range(7))
@pytest.mark.parametrize('name', ['is_polynomial', 'is_rational_function', 'is_algebraic_expr'])
def test_algebraic_queries(index, name):
    x, a, _, es = expressions()
    e = es[index]
    for variables in [(), (x,), (a,)]:
        assert getattr(e, name)(*variables) == getattr(e.to_sympy(), name)(
            *(v.to_sympy() for v in variables))


@pytest.mark.parametrize('index', range(7))
def test_coefficient_and_independence(index):
    from test_matrixexpr_protocols import convert
    x, a, b, es = expressions()
    e = es[index]
    for target in [a, b, a*b, x]:
        for right in [False, True]:
            assert convert(e.coeff(target, right=right)) == e.to_sympy().coeff(target.to_sympy(), right=right)
        assert convert(e.as_independent(target)) == e.to_sympy().as_independent(target.to_sympy())
    assert e.is_constant(x) == e.to_sympy().is_constant(x.to_sympy())


@pytest.mark.parametrize('index', range(7))
def test_series_taylor_and_integral(index):
    x, _, _, es = expressions()
    e = es[index]
    series = e.series(x, 0, 5).removeO().as_explicit()
    reference = e.to_sympy().as_explicit()
    for row in range(2):
        for col in range(2):
            expected = reference[row, col].series(x.to_sympy(), 0, 5).removeO()
            assert sp.simplify(series[row, col].to_sympy()-expected) == 0
    for n in range(4):
        term = e.taylor_term(n, x).as_explicit()
        for row in range(2):
            for col in range(2):
                expected = sp.diff(reference[row, col], x.to_sympy(), n).subs(x.to_sympy(), 0)*x.to_sympy()**n/sp.factorial(n)
                assert sp.simplify(term[row, col].to_sympy()-expected) == 0
    integral = e.integrate((x, 0, 1)).as_explicit()
    expected = reference.applyfunc(lambda v: sp.integrate(v, (x.to_sympy(), 0, 1)))
    assert (integral.to_sympy()-expected).applyfunc(sp.simplify) == sp.zeros(2)


def test_limits_and_asymptotic_series():
    x = sx.Symbol('x', positive=True)
    a = sx.MatrixSymbol('A', 2, 2)
    e = (1+1/x+1/x**2)*a
    assert e.limit(x, sx.oo) == a
    assert e.aseries(x, n=3).removeO().doit() == e
    assert (x*x*a).dir(x, 2) == 8*a


def test_function_matrix_entrywise_integration():
    i, j = sx.symbols('i j', integer=True)
    x = sx.Symbol('x')
    e = sx.FunctionMatrix(2, 3, sx.Lambda((i, j), i*x+j*x*x))
    integral = e.integrate((x, 0, 1))
    assert isinstance(integral, sx.FunctionMatrix)
    assert integral.as_explicit() == sx.Matrix(2, 3, lambda i, j: sx.Rational(i, 2)+sx.Rational(j, 3))
    with pytest.raises(ValueError):
        e.integrate(i)


def test_matrix_equality_and_scalar_function_queries():
    x, a, b, _ = expressions()
    assert ((sx.sin(x)**2+sx.cos(x)**2)*a).equals(a) is True
    assert a.equals(b) is False
    assert a.equals(sx.MatrixSymbol('C', 3, 3)) is False
    assert (sx.sin(x)*a).has(sx.sin)
    assert (sx.sin(x)*a).find(sx.sin) == {sx.sin(x)}
    assert not a.has('A')


def test_consistent_repeated_wildcard_bindings():
    _, a, b, _ = expressions()
    w, z = sx.Wild('w'), sx.Wild('z')
    pattern = w*a+w*b
    assert pattern.matches(3*a+4*b) is None
    assert (w*a+z*b).matches(3*a+4*b) == {w: sx.Integer(3), z: sx.Integer(4)}


@pytest.mark.parametrize('binder', ['Lambda', 'Integral', 'Sum', 'Product'])
def test_scalar_free_queries_respect_bound_variables(binder):
    x, y = sx.symbols('x y')
    f = sx.Function('f')
    if binder == 'Lambda':
        expression = sx.Lambda(x, f(x)+y)
    else:
        expression = getattr(sx, binder)(f(x)+y, (x, 0, 3))
    reference = expression.to_sympy()
    for needle in [x, y, f(x), sx.Integer(3)]:
        assert expression.has_free(needle) == reference.has_free(needle.to_sympy())
        assert expression.has_xfree({needle}) == reference.has_xfree({needle.to_sympy()})
