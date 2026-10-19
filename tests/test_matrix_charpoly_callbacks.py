"""Characteristic polynomial callbacks see native, unsimplified coefficients."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def transform(module,mode,x,y):
    if mode=='identity':return lambda value:value
    if mode=='simplify':return module.simplify
    if mode=='substitute':return lambda value:value.subs(x,0)
    if mode=='half':return lambda value:value/2
    if mode=='zero':return lambda value:0
    if mode=='symbol':return lambda value:y
    raise AssertionError(mode)


@pytest.mark.parametrize('kind',['Matrix','ImmutableMatrix','SparseMatrix','ImmutableSparseMatrix'])
@pytest.mark.parametrize('mode',['identity','simplify','substitute','half','zero','symbol'])
@pytest.mark.parametrize('source',[
    '[[1,2],[3,4]]',
    '[[x,1],[2,x]]',
    '[[sin(x)**2+cos(x)**2,1],[1,2]]',
])
def test_coefficient_transformations_match_reference(kind,mode,source):
    x,y,t=sx.symbols('x y t');X,Y,T=sp.symbols('x y t')
    matrix=getattr(sx,kind)(eval(source,{'x':x,'sin':sx.sin,'cos':sx.cos}))
    reference=getattr(sp,kind)(eval(source,{'x':X,'sin':sp.sin,'cos':sp.cos}))
    original=matrix.copy()
    actual=matrix.charpoly(t,simplify=transform(sx,mode,x,y))
    expected=reference.charpoly(T,simplify=transform(sp,mode,X,Y))
    assert actual.gen.to_sympy()==expected.gen
    assert len(actual.all_coeffs())==len(expected.all_coeffs())
    for coefficient,foreign in zip(actual.all_coeffs(),expected.all_coeffs()):
        assert sp.expand(coefficient.to_sympy()-foreign)==0
    assert str(actual.domain)==str(expected.domain)
    degree=actual.degree()
    assert (degree.to_sympy() if hasattr(degree,'to_sympy') else degree)==expected.degree()
    assert matrix==original


def test_callbacks_receive_unsimplified_coefficients_in_descending_order():
    x=sx.Symbol('x');expression=sx.sin(x)**2+sx.cos(x)**2
    seen=[]
    result=sx.Matrix([[expression,0],[0,2]]).charpoly(simplify=lambda value:seen.append(value) or value)
    assert seen==result.all_coeffs()
    assert any(value.has(sx.sin(x)) for value in seen)
    assert result!=sx.Matrix([[expression,0],[0,2]]).charpoly()


@pytest.mark.parametrize('dummy',[False,True])
@pytest.mark.parametrize('collision',[False,True])
def test_generator_identity_survives_callback_conversion(dummy,collision):
    x=(sx.Dummy if dummy else sx.Symbol)('x',real=True)
    y=sx.Symbol('y')
    matrix=sx.Matrix([[x if collision else y,1],[0,2]])
    result=matrix.charpoly(x,simplify=lambda value:value)
    if collision:assert result.gen!=x and x in result.as_expr().free_symbols
    else:assert result.gen==x


def test_callback_can_return_the_generator_as_a_formal_coefficient():
    x=sx.Symbol('x');result=sx.eye(2).charpoly(x,simplify=lambda value:x)
    assert result.all_coeffs()==[x,x,x] and result.degree()==2
    assert str(result.domain)=='ZZ[x]'


def test_invalid_callbacks_and_callback_exceptions():
    class Marker(Exception):pass
    failure=Marker('coefficient failure');calls=[]
    def stop(value):
        calls.append(value)
        raise failure
    matrix=sx.eye(2)
    with pytest.raises(Marker) as caught:matrix.charpoly(simplify=stop)
    assert caught.value is failure and calls==[sx.S.One]
    for callback in [True,False,42,sx.Symbol('x')]:
        with pytest.raises(TypeError):matrix.charpoly(simplify=callback)
    with pytest.raises(TypeError):matrix.charpoly(simplify=lambda value:None)
    with pytest.raises(sx.NonSquareMatrixError):sx.zeros(2,3).charpoly(simplify=lambda value:value)


def test_native_lambda_and_zero_polynomial():
    x=sx.Symbol('x')
    assert sx.eye(2).charpoly(simplify=sx.Lambda(x,x))==sx.eye(2).charpoly()
    polynomial=sx.eye(2).charpoly(simplify=lambda value:0)
    assert polynomial.all_coeffs()==[0] and polynomial.degree()==-sx.oo


def test_pure_polynomial_equality_unifies_exact_coefficient_domains():
    rational=sx.Matrix([[sx.Rational(1,2),sx.Rational(1,2)],[sx.Rational(1,2),sx.Rational(1,2)]]).charpoly()
    integer=sx.diag(0,1).charpoly()
    assert rational.all_coeffs()==integer.all_coeffs()
    reference_rational=sp.Matrix([[sp.Rational(1,2)]*2]*2).charpoly()
    reference_integer=sp.diag(0,1).charpoly()
    assert (rational==integer)==(reference_rational==reference_integer)
    assert hash(rational)==hash(integer) and len({rational,integer})==1
    x,y=sx.symbols('x y')
    assert sx.eye(2).charpoly(x)==sx.eye(2).charpoly(y)
    assert hash(sx.eye(2).charpoly(x))==hash(sx.eye(2).charpoly(y))


def test_charpoly_callbacks_do_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x=symbols('x')
for cls in [Matrix,ImmutableMatrix,SparseMatrix,ImmutableSparseMatrix]:
    matrix=cls([[x,1],[2,x]])
    polynomial=matrix.charpoly(simplify=lambda value:value.subs(x,0))
    assert polynomial.all_coeffs()==[1,0,-2]
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=30)
