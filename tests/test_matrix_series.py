"""Exact scalar series kernels with noncommutative matrix coefficients."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from symox.series.matrix import MatrixFormalPowerSeries


@pytest.mark.parametrize('center,direction', [(0, 1), (0, -1), (1, 1), (1, -1)])
@pytest.mark.parametrize('kind', ['exp', 'sin', 'rational', 'polynomial'])
def test_formal_coefficients_and_truncation(center, direction, kind):
    x = sx.Symbol('x')
    a, b = sx.MatrixSymbol('A', 2, 2), sx.MatrixSymbol('B', 2, 2)
    scalar = {'exp': sx.exp(x), 'sin': sx.sin(x), 'rational': 1/(2-x),
              'polynomial': x*x+1}[kind]
    expression = scalar*a+x*b
    formal = expression.fps(x, center, direction)
    # The matrix series keeps its own lazy class (no reference counterpart).
    assert isinstance(formal, MatrixFormalPowerSeries)
    for n in range(5):
        expected = scalar.diff(x, n).subs(x, center)*(direction**n/sx.factorial(n))*a
        if n == 0:
            expected += center*b
        elif n == 1:
            expected += direction*b
        assert formal.coefficient(n).expand().doit() == expected.expand().doit()
    assert formal.truncate(4).removeO().doit() == formal.polynomial(4).doit()
    explicit = formal.polynomial(4).as_explicit().to_sympy()
    expected = expression.as_explicit().to_sympy().applyfunc(
        lambda entry: entry.series(x.to_sympy(), center, 4).removeO())
    assert (explicit-expected).applyfunc(sp.simplify) == sp.zeros(2)
    assert pickle.loads(pickle.dumps(formal)).coefficient(2) == formal.coefficient(2)


def test_matrix_formal_product_preserves_factor_order():
    x = sx.Symbol('x')
    a, b, c, d = [sx.MatrixSymbol(name, 2, 2) for name in 'ABCD']
    left, right = (a+x*b).fps(x), (c+x*d).fps(x)
    product = left.product(right, n=4)
    assert product.coefficient(0) == a*c
    assert product.coefficient(1) == a*d+b*c
    assert product.coefficient(2) == b*d
    assert product.coefficient(3) == sx.ZeroMatrix(2, 2)
    assert left.product(right).coefficient(1) != right.product(left).coefficient(1)
    assert product.polynomial(4).expand() == ((a+x*b)*(c+x*d)).expand()


@pytest.mark.parametrize('direction', [1, -1])
def test_matrix_inverse_series_identity(direction):
    x = sx.Symbol('x')
    a, b = sx.MatrixSymbol('A', 2, 2), sx.MatrixSymbol('B', 2, 2)
    original = (a+x*b).fps(x, dir=direction)
    inverse = original.inverse()
    assert inverse.coefficient(0) == a.I
    assert inverse.coefficient(1) == -direction*a.I*b*a.I
    assert inverse.coefficient(2) == a.I*b*a.I*b*a.I
    for product in [original.product(inverse), inverse.product(original)]:
        assert product.coefficient(0).doit() == sx.Identity(2)
        for n in range(1, 5):
            assert product.coefficient(n).expand().doit() == sx.ZeroMatrix(2, 2)


def test_scalar_inner_composition_and_formal_options():
    x = sx.Symbol('x')
    n = sx.Symbol('n', integer=True, positive=True)
    a = sx.MatrixSymbol('A', n, n)
    formal = (sx.exp(x)*a).fps()
    assert formal.x == x
    inner = sx.fps(sx.sin(x), x)
    composed = formal.compose(inner)
    assert composed.polynomial(4).expand() == (a+ x*a+x*x/2*a).expand()
    assert a.fps() == a
    assert sx.fps(sx.exp(x)*a, x).coefficient(2) == a/2
    assert (sx.sqrt(x)*a).fps(x) == sx.sqrt(x)*a
    assert ((1/x)*a).fps(x) == (1/x)*a
    assert (sx.exp(x)*a).fps(x, hyper=False, rational=False) == sx.exp(x)*a
    with pytest.raises(ValueError):
        formal.compose((x*a).fps(x))
    with pytest.raises(ValueError):
        (x*a).fps(x).inverse()


def test_formal_differentiation_and_integration():
    x = sx.Symbol('x')
    a = sx.MatrixSymbol('A', 2, 2)
    formal = (sx.exp(x)*a).fps(x)
    assert formal.diff().coefficient(3) == a/6
    assert formal.integrate().coefficient(3) == a/6


@pytest.mark.parametrize('kind', ['x', 'square', 'trig', 'shifted_trig'])
def test_matrix_fourier_coefficients(kind):
    x = sx.Symbol('x')
    a, b = sx.MatrixSymbol('A', 2, 2), sx.MatrixSymbol('B', 2, 2)
    scalar = {'x': x, 'square': x*x, 'trig': 2+sx.sin(x)+3*sx.cos(2*x),
              'shifted_trig': sx.sin(2*x+sx.pi/3)}[kind]
    matrix = scalar*a+sx.cos(x)*b
    series = matrix.fourier_series((x, -sx.pi, sx.pi))
    assert isinstance(series, sx.FourierSeries)
    reference_a = sp.fourier_series(scalar.to_sympy(), (x.to_sympy(), -sp.pi, sp.pi))
    reference_b = sp.fourier_series(sp.cos(x.to_sympy()), (x.to_sympy(), -sp.pi, sp.pi))
    assert series.a0.to_sympy() == reference_a.a0*a.to_sympy()+reference_b.a0*b.to_sympy()
    for harmonic in range(1, 5):
        for sine in [False, True]:
            wave = sx.sin if sine else sx.cos
            wave_ref = sp.sin if sine else sp.cos
            coefficient = sp.integrate(scalar.to_sympy()*wave_ref(harmonic*x.to_sympy()),
                                       (x.to_sympy(), -sp.pi, sp.pi))/sp.pi
            expected = coefficient*a.to_sympy()
            if not sine and harmonic == 1:
                expected += b.to_sympy()
            actual = series._coefficient(harmonic, sine)
            assert sp.simplify(actual.to_sympy()-expected) == sp.ZeroMatrix(2, 2)


def test_finite_matrix_fourier_truncation_and_transforms():
    x = sx.Symbol('x')
    a, b, c = [sx.MatrixSymbol(name, 2, 2) for name in 'ABC']
    matrix = a+sx.sin(x)*b+sx.cos(2*x)*c
    series = sx.fourier_series(matrix, (x, -sx.pi, sx.pi))
    assert series.truncate(20) == matrix
    assert list(series) == [a, sx.sin(x)*b, sx.cos(2*x)*c]
    assert series.scale(2).truncate(20) == 2*matrix
    assert series.shift(a).a0 == 2*a
    assert series.shiftx(1).truncate(20) == matrix.subs(x, x+1)
    assert series.scalex(2).truncate(20) == matrix.subs(x, 2*x)
    sigma = series.sigma_approximation(3)
    assert sigma == a+sx.sinc(sx.pi/3)*sx.sin(x)*b+sx.sinc(2*sx.pi/3)*sx.cos(2*x)*c
    assert pickle.loads(pickle.dumps(series)).truncate(20) == matrix
    with pytest.raises(ValueError):
        series.shift(sx.MatrixSymbol('D', 3, 3))
    with pytest.raises(ValueError):
        matrix.fourier_series((x, 0, 0))


def test_matrix_series_without_sympy():
    code = '''
import sys
sys.modules['sympy']=None
import symox as s
x=s.Symbol('x'); a=s.MatrixSymbol('A',2,2); b=s.MatrixSymbol('B',2,2)
f=(a+x*b).fps(x)
assert f.coefficient(1)==b
assert f.inverse().coefficient(1)==-a.I*b*a.I
assert (a+s.sin(x)*b).fourier_series((x,-s.pi,s.pi)).truncate(10)==a+s.sin(x)*b
'''
    subprocess.run([sys.executable, '-c', code], check=True, timeout=30)
