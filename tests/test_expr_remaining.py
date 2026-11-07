"""Behavioral checks for every formerly absent public Expr method."""
import itertools
import pickle
import subprocess
import sys

import pytest
import sympy as sp
from sympy.series.formal import FormalPowerSeries
import symox as sx


def equivalent(actual, expected):
    back = actual.to_sympy()
    assert back == expected or sp.simplify(back-expected) == 0, (actual, expected)


@pytest.mark.parametrize('source', [
    'x+y', 'Integral(x, (x,0,y))', 'Integral(x, (x,0,x))',
    'Integral(r, (r,x))', 'Integral(x, x)',
    'Integral(x*y, (x,0,y), (y,0,x))',
    'Integral(Integral(x*y,(x,0,y)),(y,0,1))',
    'Sum(x*y,(x,0,y))', 'Product(x+y,(x,1,y))',
    'Lambda(x,x+y)', 'Lambda((x,y),x+y)',
    'Lambda(x,x+Symbol("_0"))',
    'Integral(x,(x,0,1))+Integral(y,(y,0,1))',
    'Limit(sin(x),x,0)', 'Subs(x+y,x,0)',
])
def test_as_dummy_and_canonical_variables(source):
    ours = dict(vars(sx), x=sx.Symbol('x'), y=sx.Symbol('y'), r=sx.Symbol('r', real=True))
    theirs = dict(vars(sp), x=sp.Symbol('x'), y=sp.Symbol('y'), r=sp.Symbol('r', real=True))
    actual, expected = eval(source, ours), eval(source, theirs)
    converted = {k.to_sympy(): v.to_sympy() for k, v in actual.canonical_variables.items()}
    assert converted == expected.canonical_variables
    normalized = actual.as_dummy()
    assert normalized.to_sympy() == expected.as_dummy()
    assert normalized.free_symbols == actual.free_symbols
    assert normalized.as_dummy() == normalized


def test_symbol_as_dummy_drops_facts_and_allocates_identity():
    x = sx.Symbol('x', real=True, positive=True)
    first, second = x.as_dummy(), x.as_dummy()
    assert first.is_Dummy and second.is_Dummy and first != second
    assert first.name == 'x' and first.is_real is None and first.is_positive is None
    assert first.assumptions0 == {'commutative': True}
    assert pickle.loads(pickle.dumps(first)) == first


@pytest.mark.parametrize('constructor,values', [
    ('Expr', []), ('Expr', [1]), ('Expr', [1,2,3]),
    ('Add', []), ('Add', [1,2,3]), ('Mul', []), ('Mul', [2,3,4]),
    ('Pow', [2,3]), ('Tuple', []), ('Tuple', [1,2,3]),
    ('Symbol', ['x']),
])
def test_fromiter_class_constructor(constructor, values):
    actual = getattr(sx, constructor).fromiter(iter(values))
    # Native DAG leaves normalize Python numeric inputs into symbolic numbers.
    reference_values = list(map(sp.sympify, values)) if constructor == 'Expr' else values
    expected = getattr(sp, constructor).fromiter(iter(reference_values))
    assert actual.to_sympy() == expected
    if constructor == 'Expr':
        assert actual.is_Function is False
    assert pickle.loads(pickle.dumps(actual)) == actual


@pytest.mark.parametrize('source,values', [
    ('x+y', [1,2]), ('x*y', [2,3]), ('x**2', [2,3]),
    ('sin(x)', [0]), ('Symbol("x")', ['y']), ('Tuple(x,y)', [1,2,3]),
])
def test_fromiter_instance_constructor(source, values):
    ours = dict(vars(sx), x=sx.Symbol('x'), y=sx.Symbol('y'))
    theirs = dict(vars(sp), x=sp.Symbol('x'), y=sp.Symbol('y'))
    actual = eval(source, ours).fromiter(iter(values))
    expected = eval(source, theirs).fromiter(iter(values))
    assert actual.to_sympy() == expected


def test_generic_fromiter_keeps_all_symbolic_dependencies_in_derivatives():
    x,y=sx.symbols('x y')
    value=sx.Expr.fromiter([x,y])
    X,Y=sp.symbols('x y')
    expected=sp.Expr.fromiter([X,Y])
    assert value.is_number is False
    for variable,reference in [(x,X),(y,Y)]:
        assert value.diff(variable).to_sympy() == expected.diff(reference)


@pytest.mark.parametrize('source', [
    'exp_polar(I*pi)', 'exp_polar(2*I*pi)', 'exp_polar(-I*pi)',
    'exp_polar(3*I*pi+x)', 'y*exp_polar(-5*I*pi)*exp_polar(3*I*pi+2*pi*x)',
    'exp_polar(-I*pi/2)', 'exp_polar(0)',
    'exp_polar(2*I*pi)**2', 'exp_polar(3*I*pi)**Rational(1,2)',
    'x+1', 'exp(2*I*pi)',
])
@pytest.mark.parametrize('half', [False, True])
def test_extract_branch_factor(source, half):
    ours = dict(vars(sx), x=sx.Symbol('x'), y=sx.Symbol('y'))
    theirs = dict(vars(sp), x=sp.Symbol('x'), y=sp.Symbol('y'))
    actual, expected = eval(source, ours), eval(source, theirs)
    a, n = actual.extract_branch_factor(allow_half=half)
    b, m = expected.extract_branch_factor(allow_half=half)
    assert a.to_sympy() == b
    assert n.to_sympy() == m
    assert pickle.loads(pickle.dumps(actual)) == actual


def test_polar_node_identity_derivative_and_interchange():
    x = sx.Symbol('x')
    polar = sx.exp_polar(x)
    assert polar != sx.Function('exp_polar')(x)
    assert polar.is_polar is True
    assert polar.diff(x) == polar
    assert polar.func(*polar.args) == polar
    assert sx.from_sympy(polar.to_sympy()) == polar
    assert sx.sympify(sx.srepr(polar)) == polar
    with pytest.raises((ValueError,TypeError,NotImplementedError)):
        sx.lambdify(x,polar,'numpy')


@pytest.mark.parametrize('source', [
    '1/x', '1/(1+x)', 'sin(1/x)', 'exp(1/x)', 'log(x)', 'sqrt(x)',
    'exp(-x)', 'exp(-x)+1/x', 'log(x)+exp(-x)',
    'sin(1/x+exp(-x))-sin(1/x)',
])
@pytest.mark.parametrize('hir', [False, True])
def test_asymptotic_series(source, hir):
    x, X = sx.Symbol('x'), sp.Symbol('x')
    actual = eval(source, dict(vars(sx), x=x)).aseries(x, n=4, hir=hir)
    expected = eval(source, dict(vars(sp), x=X)).aseries(X, n=4, hir=hir)
    assert actual.to_sympy() == expected


@pytest.mark.parametrize('point', [0, 1, -1, sp.oo, sp.zoo])
@pytest.mark.parametrize('source', [
    'x', 'y', '0', 'oo', 'zoo', 'nan', 'x**2+1/x',
    'sqrt(x)', 'x**log(3)', 'exp(x)', 'exp(1/x)',
    'sin(x)', 'sin(1/x)*x**2', 'cos(x)', 'tan(x)', 'log(x)',
    'log(x)+log(x+1)', 'sin(log(x))', 'gamma(x)', 'f(x)', 'f(y)',
    'asin(x)', 'acos(x)', 'atan(x)', 'acot(x)',
    'sinh(x)', 'cosh(x)', 'tanh(x)', 'asinh(x)', 'acosh(x)', 'atanh(x)',
])
def test_local_meromorphic(source, point):
    x, X = sx.Symbol('x'), sp.Symbol('x')
    ours = dict(vars(sx), x=x, y=sx.Symbol('y'), f=sx.Function('f'))
    theirs = dict(vars(sp), x=X, y=sp.Symbol('y'), f=sp.Function('f'))
    actual, expected = sx.sympify(eval(source, ours)), sp.sympify(eval(source, theirs))
    assert actual.is_meromorphic(x, sx.from_sympy(point)) is expected.is_meromorphic(X, point)


@pytest.mark.parametrize('source', ['sin(x)', 'cos(x)', 'exp(x)', '1/(1-x)',
                                   'log(1+x)', 'atan(x)', 'x**2+1'])
@pytest.mark.parametrize('point,direction', [(0,1),(0,-1),(1,1)])
def test_formal_power_series(source, point, direction):
    # Same series (and same coefficients) as the reference fps, including
    # its nonzero-centre forms; a held input stays held.
    x, X = sx.Symbol('x'), sp.Symbol('x')
    function = eval(source, dict(vars(sx), x=x))
    reference = eval(source, dict(vars(sp), x=X))
    native = function.fps(x, point, direction)
    expected = reference.fps(X, point, direction)
    if source == 'x**2+1' and point == 1:
        # Deliberate deviation: the reference reads the coefficients of the
        # unexpanded shifted polynomial (x + 1)**2 + 1 and returns 1; the
        # native series keeps the Taylor coefficients.
        assert native.truncate(7).to_sympy() == \
            2 + 2*(X - 1) + (X - 1)**2 + sp.Order((X - 1)**7, (X, 1))
        assert [native[i].to_sympy() for i in range(5)] == [2, 2*(X - 1), (X - 1)**2, 0, 0]
        return
    if not isinstance(expected, FormalPowerSeries):
        assert native.to_sympy() == expected
        return
    assert isinstance(native, sx.FormalPowerSeries)
    assert native.x == x and native.x0 == point and native.dir == direction
    assert native.truncate(7).to_sympy() == expected.truncate(7)
    assert native.polynomial(7).to_sympy() == expected.polynomial(7)
    for index in range(5):
        assert native[index].to_sympy() == expected[index]
    assert [t.to_sympy() for t in itertools.islice(native.truncate(None), 4)] == \
        list(itertools.islice(expected.truncate(None), 4))
    assert pickle.loads(pickle.dumps(native)).truncate(7) == native.truncate(7)


def test_formal_series_composition_product_and_inverse():
    x, X = sx.Symbol('x'), sp.Symbol('x')
    a, b = sx.fps(sx.exp(x), x), sx.fps(sx.sin(x), x)
    A, B = sp.fps(sp.exp(X), X), sp.fps(sp.sin(X), X)
    # The reference product convolves by a symbolic FFT (unsimplified
    # radicals); the native one is exact: the same values.
    difference = sp.expand(a.product(b, n=7).truncate().removeO().to_sympy()
                           - A.product(B, n=7).truncate().removeO())
    assert all(abs(sp.N(difference.coeff(X, i), 50)) < 1e-40 for i in range(7))
    assert a.compose(b, n=7).truncate().to_sympy() == A.compose(B, n=7).truncate()
    assert a.inverse(n=7).truncate().to_sympy() == A.inverse(n=7).truncate()
    assert a.diff().truncate().to_sympy() == A.diff(X).truncate()
    assert b.integrate().truncate().to_sympy() == B.integrate(X).truncate()
    for source in [1/x, sx.sqrt(x), sx.log(x), sx.exp(1/x)]:
        assert sx.fps(source, x) == source


def test_asymptotic_normalization_preserves_the_exact_exponent():
    x = sx.Symbol('x')
    f = sx.exp(sx.exp(x)/(1-1/x))
    normalized = f.aseries(x, bound=3)
    assert normalized.is_Mul
    exponent = sx.Add(*(factor.args[0] for factor in normalized.args))
    assert sp.simplify(exponent.to_sympy()-f.args[0].to_sympy()) == 0


@pytest.mark.parametrize('source', ['x', 'x**2', 'x**3', 'Abs(x)', 'sign(x)',
                                   'sin(x)', '3*cos(2*x)+sin(3*x)+2',
                                   'Piecewise((0,x<0),(1,True))'])
@pytest.mark.parametrize('interval', ['symmetric', 'positive'])
def test_fourier_series_coefficients_truncation_and_sigma(source, interval):
    x, X = sx.Symbol('x'), sp.Symbol('x')
    limits = (x, -sx.pi, sx.pi) if interval == 'symmetric' else (x, 0, 2*sx.pi)
    reference_limits = (X, -sp.pi, sp.pi) if interval == 'symmetric' else (X, 0, 2*sp.pi)
    actual = eval(source, dict(vars(sx), x=x)).fourier_series(limits)
    expected = sp.fourier_series(eval(source, dict(vars(sp), x=X)), reference_limits)
    assert isinstance(actual, sx.FourierSeries)
    equivalent(actual.a0, expected.a0)
    for i in range(4):
        try:
            term = expected[i]
        except IndexError:
            with pytest.raises(IndexError):
                actual[i]
        else:
            equivalent(actual[i], term)
    constant_interval = interval == 'positive' and source in ('sign(x)', 'Piecewise((0,x<0),(1,True))')
    # SymPy's nonzero-term iterator does not terminate for these functions,
    # which equal one throughout the open interval. Verify the exact constant.
    equivalent(actual.truncate(3), sp.S.One if constant_interval else expected.truncate(3))
    # SymPy's finite-series sigma helper raises for n above its index bound;
    # evaluate the defining finite sum with zero missing harmonics instead.
    sigma = expected.a0
    for k in range(1,4):
        try:
            term = expected[k]
        except IndexError:
            term = sp.S.Zero
        sigma += term*sp.sinc(sp.pi*k/4)
    equivalent(actual.sigma_approximation(4), sigma)
    for method, value in [('shift',1),('shiftx',1),('scale',2),('scalex',2)]:
        changed = getattr(actual, method)(value)
        if constant_interval:
            reference_value = sp.Integer(2) if method in ('shift','scale') else sp.S.One
        else:
            reference_value = getattr(expected, method)(value).truncate(3)
        equivalent(changed.truncate(3), reference_value)
    assert pickle.loads(pickle.dumps(actual)).truncate(3) == actual.truncate(3)


def test_errors_and_constant_series():
    x = sx.Symbol('x')
    with pytest.raises(TypeError):
        sx.sin(x).is_meromorphic(x+1, 0)
    with pytest.raises(ValueError):
        sx.sin(x).fps(x, dir='invalid')
    with pytest.raises(NotImplementedError):  # reference: multivariate
        sx.fps(x+sx.Symbol('y'))
    with pytest.raises(ValueError):
        sx.fourier_series(x, (x,0,0))
    assert sx.Integer(3).fps(x) == 3
    assert sx.Integer(3).fourier_series((x,-sx.pi,sx.pi)) == 3


def test_undefined_heads_named_like_builtins_keep_their_identity():
    x = sx.Symbol('x')
    wave = sx.Function('sin')(x).fourier_series((x,-sx.pi,sx.pi))
    assert not isinstance(wave,sx.FiniteFourierSeries)
    with pytest.raises(NotImplementedError):
        wave.truncate(1)
    with pytest.raises((NotImplementedError,sx.PoleError)):
        sx.Function('exp')(-x).aseries(x)


def test_fourier_does_not_guess_symbolic_or_complex_branch_signs():
    x=sx.Symbol('x')
    a=sx.Symbol('a',positive=True)
    for function in (sx.Abs(x+a),sx.Abs(1+sx.I*x)):
        series=function.fourier_series((x,-sx.pi,sx.pi))
        with pytest.raises(NotImplementedError):
            series.truncate(1)


def test_all_remaining_methods_without_sympy_runtime():
    subprocess.run([sys.executable, '-c', '''
import sys
sys.modules['sympy'] = None
import symox as s
x=s.Symbol('x')
assert s.Integral(x,(x,0,x)).as_dummy().free_symbols == {x}
assert s.Add.fromiter([1,2,3]) == 6
assert s.Expr.fromiter([1,2]).args == (1,2)
assert s.exp_polar(2*s.pi*s.I).extract_branch_factor() == (s.Integer(1),s.Integer(1))
assert s.log(x).is_meromorphic(x,0) is False
assert (1/x).aseries(x) == 1/x
assert s.sin(x).fps(x).truncate(4) == x-x**3/6+s.O(x**4,x)
assert x.fourier_series((x,-s.pi,s.pi)).truncate(2) == 2*s.sin(x)-s.sin(2*x)
'''], check=True, timeout=30)
