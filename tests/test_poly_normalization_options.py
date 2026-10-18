"""Check polynomial option contracts as well as exact domains and values."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def polynomial_pair(case):
    x, y, a, b = sx.symbols('x y a b')
    X, Y, A, B = sp.symbols('x y a b')
    values = [
        (2*x+4, (x,), sx.ZZ, 2*X+4, (X,), sp.ZZ),
        (2*x+3, (x,), sx.ZZ, 2*X+3, (X,), sp.ZZ),
        (0, (x,), sx.ZZ, 0, (X,), sp.ZZ),
        (x/2+sx.Rational(1, 3), (x,), sx.QQ, X/2+sp.Rational(1, 3), (X,), sp.QQ),
        (0, (x,), sx.QQ, 0, (X,), sp.QQ),
        (2*x+4, (x,), sx.GF(5), 2*X+4, (X,), sp.GF(5)),
        (a*x+a, (x,), sx.ZZ.poly_ring(a), A*X+A, (X,), sp.ZZ.poly_ring(A)),
        (a*x+1, (x,), sx.ZZ.poly_ring(a), A*X+1, (X,), sp.ZZ.poly_ring(A)),
        (2*a*x+a, (x,), sx.QQ.poly_ring(a), 2*A*X+A, (X,), sp.QQ.poly_ring(A)),
        (x/(2*a+2)+sx.Rational(1, 3), (x,), sx.ZZ.frac_field(a),
         X/(2*A+2)+sp.Rational(1, 3), (X,), sp.ZZ.frac_field(A)),
        (x/(2*a+2)+sx.Rational(1, 3), (x,), sx.QQ.frac_field(a),
         X/(2*A+2)+sp.Rational(1, 3), (X,), sp.QQ.frac_field(A)),
        (x/2+1, (x,), sx.QQ.poly_ring(a), X/2+1, (X,), sp.QQ.poly_ring(A)),
        (x/2+1, (x,), sx.QQ.frac_field(a), X/2+1, (X,), sp.QQ.frac_field(A)),
        (x/(2*a**2+3*b**3)+1/(a+1), (x,), sx.QQ.frac_field(a,b,order='grlex'),
         X/(2*A**2+3*B**3)+1/(A+1), (X,), sp.QQ.frac_field(A,B,order='grlex')),
        ((sx.sqrt(2)+1)*x/3+1, (x,), sx.QQ.algebraic_field(sx.sqrt(2)),
         (sp.sqrt(2)+1)*X/3+1, (X,), sp.QQ.algebraic_field(sp.sqrt(2))),
        (2*x*y+4*x+6*y, (x,y), sx.ZZ, 2*X*Y+4*X+6*Y, (X,Y), sp.ZZ),
    ]
    expression, gens, domain, expected, sp_gens, sp_domain = values[case]
    return sx.Poly(expression, *gens, domain=domain), sp.Poly(expected, *sp_gens, domain=sp_domain)


def assert_same_polynomial(actual, expected):
    assert str(actual.domain) == str(expected.domain)
    assert sp.cancel(actual.as_expr().to_sympy()-expected.as_expr()) == 0
    if actual.domain.is_Composite:
        assert actual.domain.order == str(expected.domain.order)


@pytest.mark.parametrize('case', range(16))
@pytest.mark.parametrize('auto', [False, True])
def test_monic_option_and_zero_domain(case, auto):
    polynomial, expected = polynomial_pair(case)
    try:
        expected_result = expected.monic(auto=auto)
    except sp.polys.polyerrors.ExactQuotientFailed:
        with pytest.raises(sx.ExactQuotientFailed):
            polynomial.monic(auto=auto)
    else:
        result = polynomial.monic(auto=auto)
        assert_same_polynomial(result, expected_result)
        if not result.is_zero:
            assert result.LC() == 1


@pytest.mark.parametrize('case', range(16))
@pytest.mark.parametrize('convert', [False, True])
def test_clear_denominators_option_and_ground_ring(case, convert):
    polynomial, expected = polynomial_pair(case)
    denominator, numerator = polynomial.clear_denoms(convert=convert)
    expected_denominator, expected_numerator = expected.clear_denoms(convert=convert)
    assert sp.cancel(denominator.to_sympy()-expected_denominator) == 0
    assert_same_polynomial(numerator, expected_numerator)
    assert sp.cancel((denominator*polynomial.as_expr()-numerator.as_expr()).to_sympy()) == 0


@pytest.mark.parametrize('case', range(8))
@pytest.mark.parametrize('domain_kind', ['ZZ','QQ','GF','ring','field','algebraic'])
def test_derivative_orders_and_signed_indices(case, domain_kind):
    x, y, a = sx.symbols('x y a')
    X, Y, A = sp.symbols('x y a')
    domains = dict(ZZ=(sx.ZZ,sp.ZZ), QQ=(sx.QQ,sp.QQ), GF=(sx.GF(5),sp.GF(5)),
        ring=(sx.ZZ.poly_ring(a),sp.ZZ.poly_ring(A)),
        field=(sx.QQ.frac_field(a),sp.QQ.frac_field(A)),
        algebraic=(sx.QQ.algebraic_field(sx.sqrt(2)),sp.QQ.algebraic_field(sp.sqrt(2))))
    domain, expected_domain = domains[domain_kind]
    polynomial = sx.Poly(3*x**4*y**3+2*x*y+1, x,y,domain=domain)
    expected = sp.Poly(3*X**4*Y**3+2*X*Y+1, X,Y,domain=expected_domain)
    specs = [(),((x,2),),((x,0),),((x,-1),),(-1,),((0,3),),((x,2),(y,2)),((x,1000000),)]
    sp_specs = [(),((X,2),),((X,0),),((X,-1),),(-1,),((0,3),),((X,2),(Y,2)),((X,1000000),)]
    assert_same_polynomial(polynomial.diff(*specs[case]),expected.diff(*sp_specs[case]))
    assert polynomial.degree(-1) == expected.degree(-1)


def test_invalid_generators_and_derivative_specs_are_checked():
    x, y = sx.symbols('x y')
    polynomial = sx.Poly(x*x+y,x,y)
    for index in [-3,2,1000]:
        for method in [polynomial.diff,polynomial.degree]:
            with pytest.raises(ValueError):
                method(index)
    for spec in [(x,), (x,2,3)]:
        with pytest.raises(ValueError):
            polynomial.diff(spec)
    with pytest.raises(TypeError):
        polynomial.diff((x,1.5))
    with pytest.raises(ValueError):
        polynomial(1,2,3)


def test_global_monic_preserves_ring_with_auto_false():
    x,a = sx.symbols('x a')
    polynomial = sx.Poly(2*x+4,x,domain=sx.ZZ)
    assert sx.monic(polynomial,auto=False).domain == sx.ZZ
    assert sx.monic(a*x+a,x,auto=False) == x+1
    with pytest.raises(sx.ExactQuotientFailed):
        sx.monic(a*x+1,x,auto=False)


def test_no_sympy_import():
    script = '''
import sys
class Block:
    def find_spec(self, fullname, *args):
        if fullname == 'sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import Poly, symbols, ZZ, QQ
x,a = symbols('x a')
assert Poly(2*x+4,x).monic(auto=False).domain == ZZ
assert Poly(x**4,x).diff((x,3)).as_expr() == 24*x
K=QQ.frac_field(a)
p=Poly(x/(2*a+2)+1,x,domain=K)
d,n=p.clear_denoms(convert=True)
assert d == a+1 and n.domain == QQ.poly_ring(a)
'''
    subprocess.run([sys.executable,'-c',script],check=True)
