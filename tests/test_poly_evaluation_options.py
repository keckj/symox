"""Domain-preserving partial evaluation and expression substitutions."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def pair(kind):
    x,y,z,a,b=sx.symbols('x y z a b')
    X,Y,Z,A,B=sp.symbols('x y z a b')
    domains=[sx.ZZ,sx.QQ,sx.GF(5),sx.ZZ.poly_ring(a,b,order='grlex'),
             sx.QQ.frac_field(a,b,order='grevlex'),sx.QQ.algebraic_field(sx.sqrt(2),alias='r')]
    reference=[sp.ZZ,sp.QQ,sp.GF(5),sp.ZZ.poly_ring(A,B,order='grlex'),
               sp.QQ.frac_field(A,B,order='grevlex'),sp.QQ.algebraic_field(sp.sqrt(2),alias='r')]
    domain,ref_domain=domains[kind],reference[kind]
    coefficient=a if kind==3 else 1/a if kind==4 else sx.sqrt(2) if kind==5 else 2
    ref_coefficient=A if kind==3 else 1/A if kind==4 else sp.sqrt(2) if kind==5 else 2
    expression=coefficient*x**3*y+3*x*y*y+z+1
    expected=ref_coefficient*X**3*Y+3*X*Y*Y+Z+1
    return sx.Poly(expression,x,y,z,domain=domain),sp.Poly(expected,X,Y,Z,domain=ref_domain)


def assert_equivalent(actual,expected):
    if isinstance(expected,sp.Poly):
        assert isinstance(actual,sx.Poly)
        assert str(actual.domain)==str(expected.domain)
        assert tuple(g.to_sympy() for g in actual.gens)==expected.gens
        assert sp.cancel(actual.as_expr().to_sympy()-expected.as_expr())==0
        if actual.domain.is_Composite:
            assert actual.domain.order==str(expected.domain.order)
    else:
        assert not isinstance(actual,sx.Poly)
        assert sp.simplify(actual.to_sympy()-expected)==0


@pytest.mark.parametrize('kind',range(6))
@pytest.mark.parametrize('count',range(4))
@pytest.mark.parametrize('container',['tuple','list','dict'])
def test_evaluation_containers_and_remaining_polynomial_domains(kind,count,container):
    polynomial,expected=pair(kind)
    values=[3,7,11][:count]
    native_values=dict(zip(polynomial.gens,values)) if container=='dict' else (
        values if container=='list' else tuple(values))
    reference_values=dict(zip(expected.gens,values)) if container=='dict' else (
        values if container=='list' else tuple(values))
    assert_equivalent(polynomial.eval(native_values),expected.eval(reference_values))
    assert_equivalent(polynomial(*values),expected(*values))


@pytest.mark.parametrize('kind',range(6))
@pytest.mark.parametrize('index',[-3,-2,-1,0,1,2])
def test_generator_and_signed_index_evaluation(kind,index):
    polynomial,expected=pair(kind)
    assert_equivalent(polynomial.eval(index,17),expected.eval(index,17))
    assert_equivalent(polynomial.eval(polynomial.gens[index],17),
                      expected.eval(expected.gens[index],17))


@pytest.mark.parametrize('kind',[0,1,3,4,5])
@pytest.mark.parametrize('auto',[False,True])
def test_rational_promotion_and_exact_auto_false(kind,auto):
    polynomial,expected=pair(kind)
    try:
        reference=expected.eval(expected.gens[0],sp.Rational(1,2),auto=auto)
    except sp.polys.polyerrors.DomainError:
        with pytest.raises(sx.DomainError):
            polynomial.eval(polynomial.gens[0],sx.Rational(1,2),auto=auto)
    else:
        assert_equivalent(polynomial.eval(polynomial.gens[0],sx.Rational(1,2),auto=auto),reference)


@pytest.mark.parametrize('point',[4,5,6,10,123456789012345678901234567890])
def test_finite_field_evaluation_reduces_the_result(point):
    x=sx.Symbol('x');X=sp.Symbol('x')
    polynomial=sx.Poly(x**3+2*x+1,x,modulus=5)
    expected=sp.Poly(X**3+2*X+1,X,modulus=5)
    assert_equivalent(polynomial.eval(point),expected.eval(point))


def test_promotion_to_symbolic_ground_ring_and_collision_checks():
    x,y,z=sx.symbols('x y z');X,Y,Z=sp.symbols('x y z')
    polynomial=sx.Poly(x*x+y,x,y)
    expected=sp.Poly(X*X+Y,X,Y)
    assert_equivalent(polynomial.eval(x,z+1),expected.eval(X,Z+1))
    for point in [x,y,x+y]:
        with pytest.raises(sx.UnificationFailed):
            polynomial.eval(x,point)
    with pytest.raises(sx.DomainError):
        polynomial.eval(x,z,auto=False)
    with pytest.raises(sx.UnificationFailed):
        sx.Poly(x+y,x,y,modulus=5).eval(x,sx.Rational(1,2))
    # The explicit retention flag applies throughout a container evaluation.
    with pytest.raises(sx.DomainError):
        polynomial.eval({x:sx.Rational(1,2)},auto=False)


@pytest.mark.parametrize('point',[sx.sqrt(2),sx.I,sx.pi,sx.sin(sx.Symbol('z'))])
def test_univariate_evaluation_can_return_general_native_scalars(point):
    x=sx.Symbol('x');X=sp.Symbol('x')
    polynomial=sx.Poly(x*x+2*x+3,x,domain=sx.QQ)
    expected=sp.Poly(X*X+2*X+3,X,domain=sp.QQ)
    assert_equivalent(polynomial.eval(point),expected.eval(point.to_sympy()))


def test_new_coefficient_backends_follow_the_reference():
    x,y=sx.symbols('x y')
    polynomial=sx.Poly(x+y,x,y)
    # The reference evaluates at I into ZZ_I.
    value=polynomial.eval(x,sx.I)
    assert str(value)=="Poly(y + I, y, domain='ZZ_I')"
    # Algebraic points without an extension land in EX, like the reference.
    value=polynomial.eval(x,sx.sqrt(2))
    assert str(value)=="Poly(y + sqrt(2), y, domain='EX')"


@pytest.mark.parametrize('case',range(8))
def test_as_expr_replaces_generators_simultaneously(case):
    x,y,z=sx.symbols('x y z');X,Y,Z=sp.symbols('x y z')
    polynomial=sx.Poly(x*x+2*x*y*y-y,x,y)
    expected=sp.Poly(X*X+2*X*Y*Y-Y,X,Y)
    args=[(),({x:5},),(5,6),({x:y,y:x},),(z+1,z-1),({x:x+y},),
          (sx.Rational(1,2),sx.sqrt(2)),(0,0)][case]
    ref_args=[(),({X:5},),(5,6),({X:Y,Y:X},),(Z+1,Z-1),({X:X+Y},),
              (sp.Rational(1,2),sp.sqrt(2)),(0,0)][case]
    # SymPy 1.14 searches its already modified generator list for each key,
    # so swapping two generators through a dictionary cancels the swap.
    reference = (expected.as_expr().subs({X:Y,Y:X},simultaneous=True) if case==3
                 else expected.as_expr(*ref_args))
    assert_equivalent(polynomial.as_expr(*args),reference)


@pytest.mark.parametrize('point',['0.0','0.1','1.23456789012345','-3.25'])
def test_numerical_evaluation_uses_field_precision(point):
    x=sx.Symbol('x');X=sp.Symbol('x')
    polynomial=sx.Poly(sx.Float('1.25')*x*x+sx.Float('0.5')*x+1,x)
    expected=sp.Poly(sp.Float('1.25')*X*X+sp.Float('0.5')*X+1,X)
    result=polynomial.eval(sx.Float(point))
    reference=expected.eval(sp.Float(point))
    assert sx.RR.almosteq(result,sx.sympify(reference))


def test_as_expr_keeps_coefficient_parameters_independent_of_the_mapping():
    x,a=sx.symbols('x a');X,A=sp.symbols('x a')
    polynomial=sx.Poly(a*x+a*x*x,x,domain=sx.ZZ.poly_ring(a))
    expected=sp.Poly(A*X+A*X*X,X,domain=sp.ZZ.poly_ring(A))
    assert_equivalent(polynomial.as_expr({x:a}),expected.as_expr({X:A}))
    with pytest.raises(sx.GeneratorsError):
        polynomial.as_expr({a:2})


def test_expression_generator_evaluation_uses_polynomial_positions():
    x,t=sx.symbols('x t')
    polynomial=sx.Poly(t*t+2*t+1,t)._with_gens(x+1)
    assert polynomial.eval(3)==16
    assert polynomial.as_expr(5)==36


def test_invalid_evaluation_arity_and_generators():
    x,y,z=sx.symbols('x y z')
    polynomial=sx.Poly(x+y,x,y)
    for values in [[1,2,3],(1,2,3)]:
        with pytest.raises(ValueError):
            polynomial.eval(values)
    for generator in [-3,2,z]:
        with pytest.raises(sx.GeneratorsError):
            polynomial.eval(generator,1)
    for args in [(1,),(1,2,3),({z:1},)]:
        with pytest.raises(sx.GeneratorsError):
            polynomial.as_expr(*args)


def test_no_sympy_import():
    script='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import Poly,symbols,ZZ
x,y=symbols('x y')
p=Poly(x*x+y,x,y,domain=ZZ)
assert p.eval((3,4))==13
assert p.eval({x:3}).gens==(y,)
assert Poly(x+1,x,modulus=5).eval(4)==0
assert p.as_expr({x:y})==y*y+y
'''
    subprocess.run([sys.executable,'-c',script],check=True)
