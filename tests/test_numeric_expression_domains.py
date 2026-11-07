"""Native real/complex precision and formal expression coefficient domains."""
from fractions import Fraction
import operator
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('precision', [53, 80, 100, 160, 200, 256, 512])
@pytest.mark.parametrize('literal', ['0.1', '-0.2', '1.23456789012345678901234567890123456789',
                                    '1e-100', '1e100'])
def test_real_decimal_input_matches_reference_bits(precision, literal):
    k, ref = sx.RealField(precision), sp.RealField(precision)
    value, expected = k(literal), ref(literal)
    assert k.precision == ref.precision and k.dps == ref.dps
    assert k.to_rational(value, limit=False) == tuple(map(int, ref.to_rational(expected, limit=False)))
    assert k.to_rational(value) == tuple(map(int, ref.to_rational(expected)))
    assert value.parent() == k and isinstance(value, k.dtype)
    assert pickle.loads(pickle.dumps(value)) == value
    assert pickle.loads(pickle.dumps(k)) == k


@pytest.mark.parametrize('precision', [53, 100, 200])
def test_numerical_closed_arithmetic_and_complex_components(precision):
    real, complex_field = sx.RealField(precision), sx.ComplexField(precision)
    assert not real.is_Exact and not complex_field.is_Exact
    assert real.is_Field and complex_field.is_Field
    a, b = real('0.1'), real('0.2')
    assert real.almosteq((a+b)-b, a)
    assert real.almosteq(a*(a**-1), real.one)
    assert ((a+b)*b).parent() == real
    z = complex_field('0.1', '0.2')
    assert z.real == a and z.imag == b
    parsed = complex_field('0.1+0.2*I')
    assert z == parsed
    assert complex_field.almosteq(z/z, complex_field.one)
    assert z.real.parent() == real and abs(z).parent() == real
    assert real.almosteq(abs(z)*abs(z), real('0.05'))
    assert complex_field.exsqrt(complex_field(-1)) == complex_field(sx.I)
    assert real.exsqrt(real(-1)) is None
    root = complex_field.exsqrt(complex_field(3, 4))
    assert complex_field.almosteq(root, complex_field(2, 1))
    assert real.get_exact() is sx.QQ
    assert complex_field.get_exact() is sx.QQ_I
    with pytest.raises(sx.CoercionFailed):
        real.convert(z)
    assert sx.QQ.convert(a) == sx.QQ(1, 10)
    assert sx.QQ_I.convert(z) == sx.QQ_I(sx.Rational(1, 10), sx.Rational(1, 5))
    with pytest.raises(sx.CoercionFailed):
        sx.QQ.convert(z)


def test_numerical_equality_and_hash_contract():
    for k in [sx.RR, sx.RealField(200), sx.CC, sx.ComplexField(200)]:
        for value in [0, 1, -2, Fraction(1, 2)]:
            a = k(value)
            assert a == value and value == a
            assert a == sx.QQ(value) and sx.QQ(value) == a
            assert hash(a) == hash(value)
            assert {a, value} == {value}
        assert k(1) == sx.RR(1)
    assert sx.RR('0.1') == 0.1 and hash(sx.RR('0.1')) == hash(0.1)
    assert sx.RealField(200)('0.1') != 0.1
    assert sx.RR('0.1') != Fraction(1, 10)
    assert sx.CC(1, 2) == complex(1, 2)
    assert hash(sx.CC(1, 2)) == hash(complex(1, 2))
    assert sx.Integer(1) == sx.RR(1) and sx.RR(1) == sx.Integer(1)
    assert sx.RR('1.5') > 1 and sx.RR('1.5') < 2
    with pytest.raises(TypeError):
        sx.CC(1, 2) < 0
    with pytest.raises(TypeError):
        operator.index(sx.RR(2))
    assert int(sx.RR('1.9')) == 1
    assert complex(sx.CC(1, 2)) == complex(1, 2)


@pytest.mark.parametrize('constructor', [sx.RealField, sx.ComplexField])
def test_numerical_constructor_options(constructor):
    assert constructor(dps=50).precision == sp.RealField(dps=50).precision
    for precision in [0, -1]:
        with pytest.raises(ValueError):
            constructor(precision)
    for precision in [1,40001]:
        assert constructor(precision).precision==precision
    with pytest.raises(TypeError):
        constructor(53, dps=20)
    with pytest.raises(TypeError):
        constructor(53.0)
    with pytest.raises(NotImplementedError):
        constructor(tol='1e-10')
    with pytest.raises(sx.CoercionFailed):
        constructor()(sx.Symbol('x'))
    with pytest.raises(sx.CoercionFailed):
        constructor()(sx.oo)


@pytest.mark.parametrize('raw', [False, True])
def test_expression_domain_arithmetic_and_structure(raw):
    k = sx.EXRAW if raw else sx.EX
    x, y = sx.symbols('x y')
    original = sx.sin(x)**2+sx.cos(x)**2
    value = k(original)
    assert value.parent() is not None and value.parent() == k
    assert (value.ex-original).trigsimp() == 0
    # SymPy 1.14 keeps the converted expression as it is (EX and EXRAW);
    # only the arithmetic normalizes.
    assert value.ex == original
    a, b = k(sx.sin(x)+y), k(x+1)
    assert a-a == k.zero and a/a == k.one
    assert (a*b)/b == a
    assert k.get_exact() is k
    assert pickle.loads(pickle.dumps(a)) == a
    assert sx.QQ.from_ExpressionRawDomain(sx.EXRAW(sx.Rational(2, 3)), sx.EXRAW) == sx.QQ(2, 3)
    assert sx.QQ.from_ExpressionDomain(sx.EX(sx.sqrt(2)), sx.EX) is None
    root = k.exsqrt(k(x))
    assert root*root == k(x)
    assert (k.to_sympy(k.from_FractionField(sx.QQ.frac_field(x)(1/x), sx.QQ.frac_field(x)))-1/x).cancel() == 0


def test_expression_sign_conventions_and_numerical_unification():
    x = sx.Symbol('x')
    for expression in [x, -x, 2*x+1, -2*x-1]:
        value = sx.EX(expression)
        reference = sp.EX.from_sympy(expression.to_sympy())
        assert sx.EX.is_positive(value) == sp.EX.is_positive(reference)
        assert sx.EX.is_negative(value) == sp.EX.is_negative(reference)
    assert sx.RR.unify(sx.QQ) == sx.RR
    assert sx.RealField(100).unify(sx.RR) == sx.RealField(100)
    assert sx.RealField(200).unify(sx.CC) == sx.ComplexField(200)
    assert sx.EX.unify(sx.RR) == sx.EX
    assert sx.EXRAW.unify(sx.EX) == sx.EXRAW
    with pytest.raises(sx.DomainError):
        sx.RR.unify(sx.GF(5))
    assert str(sx.Poly(x+1, x, domain=sx.EX)) == str(sp.Poly(sp.Symbol('x')+1, sp.Symbol('x'), domain=sp.EX))


def test_domain_constructor_facades():
    x, y = sx.symbols('x y')
    pairs = [(sx.IntegerRing(), sx.IntegerRing, sx.ZZ),
             (sx.RationalField(), sx.RationalField, sx.QQ),
             (sx.FiniteField(7), sx.FiniteField, sx.GF(7)),
             (sx.RealField(), sx.RealField, sx.RR),
             (sx.ComplexField(), sx.ComplexField, sx.CC),
             (sx.ExpressionDomain(), sx.ExpressionDomain, sx.EX),
             (sx.ExpressionRawDomain(), sx.ExpressionRawDomain, sx.EXRAW),
             (sx.GaussianIntegerRing(), sx.GaussianIntegerRing, sx.ZZ_I),
             (sx.GaussianRationalField(), sx.GaussianRationalField, sx.QQ_I),
             (sx.PolynomialRing(sx.QQ, 'x,y'), sx.PolynomialRing, sx.QQ.poly_ring(x, y)),
             (sx.FractionField(sx.ZZ, (x, y)), sx.FractionField, sx.ZZ.frac_field(x, y)),
             (sx.AlgebraicField(sx.QQ, sx.sqrt(2)), sx.AlgebraicField, sx.QQ.algebraic_field(sx.sqrt(2)))]
    for value, kind, expected in pairs:
        assert value == expected and isinstance(value, kind)
    assert sx.FF(7) == sx.GF(7)


def test_numerical_and_expression_domains_without_sympy():
    code = '''
import sys
sys.modules['sympy']=None
import symox as s
x=s.Symbol('x'); r=s.RealField(200); c=s.ComplexField(200)
assert r.to_rational(r('0.1'))==(1,10)
assert c.exsqrt(c(-1))==c(s.I)
assert c.almosteq(c('0.1','0.2')**2,c('-0.03','0.04'))
assert s.EX(s.sin(x)**2+s.cos(x)**2)!=s.EX.one
assert s.EX(s.sin(x)**2+s.cos(x)**2)-s.EX(s.sin(x)**2)==s.EX(s.cos(x)**2)
assert s.QQ.from_RealField(r('0.5'),r)==s.QQ(1,2)
'''
    subprocess.run([sys.executable, '-c', code], check=True, timeout=30)
