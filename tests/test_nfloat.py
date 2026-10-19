"""Controlled floating conversion retains exponents and function boundaries."""
from collections import defaultdict, namedtuple
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def expression(module,source):
    x,y=module.symbols('x y');f=module.Function('f')
    return module.sympify(eval(source,vars(module)|{'x':x,'y':y,'f':f}))


@pytest.mark.parametrize('digits',[5,15,40])
@pytest.mark.parametrize('exponents',[False,True])
@pytest.mark.parametrize('source',[
    '0','1','Rational(1,3)','pi','I','sqrt(2)','Rational(1,3)+I/7',
    'x**4+x/2+cos(pi/3)+1+sqrt(y)','sqrt(x)','x**Rational(2,3)',
    'x**(y/2)','x**(1+x/2)','2**x','exp(x/2)','sin(x/2)',
    'sin(cos(x/2)+Rational(1,3))','f(Rational(1,2))',
    'f(x/2)+x/3','Lambda(x,x/2)',
    'Derivative(sin(x/2),x,evaluate=False)',
    'Integral(x/2,(x,0,1))','Integral(x/2,x)',
    'Piecewise((x/2,x>0),(Rational(1,3),True))',
    'Eq(x,Rational(1,2))','Order(x*x,x)',
])
def test_floating_conversion_matches_reference(digits,exponents,source):
    actual=sx.nfloat(expression(sx,source),digits,exponent=exponents)
    expected=sp.nfloat(expression(sp,source),digits,exponent=exponents)
    assert actual.to_sympy()==expected


def test_dictionary_keys_container_types_and_default_factories():
    x=sx.Symbol('x');original={sx.Rational(1,2):[x/3,(1,sx.pi)]}
    converted=sx.nfloat(original,25)
    assert list(converted)==[sx.Rational(1,2)]
    assert converted[sx.Rational(1,2)][0]==sx.Float(sx.Rational(1,3),25)*x
    assert isinstance(converted[sx.Rational(1,2)][1],tuple)
    assert next(iter(sx.nfloat(original,dkeys=True))).is_Float
    mapping=defaultdict(int,{x:sx.Rational(1,3)})
    result=sx.nfloat(mapping)
    assert result.default_factory is int and result[sx.Symbol('missing')]==0
    Pair=namedtuple('Pair','first second')
    assert type(sx.nfloat(Pair(1,2))) is Pair
    assert isinstance(sx.nfloat({sx.Rational(1,2),sx.Rational(1,3)}),set)


@pytest.mark.parametrize('constructor',[sx.Matrix,sx.ImmutableMatrix,sx.SparseMatrix,sx.ImmutableSparseMatrix])
def test_concrete_matrices_convert_entries(constructor):
    x=sx.Symbol('x');original=constructor([[sx.Rational(1,3),sx.sqrt(x)],[1,0]])
    result=sx.nfloat(original,20)
    assert result[0,0]==sx.Float(sx.Rational(1,3),20)
    assert result[0,1]==sx.sqrt(x)
    assert result.shape==original.shape
    assert original[0,0]==sx.Rational(1,3)


@pytest.mark.parametrize('constructor',[sx.Array,sx.MutableDenseNDimArray,sx.ImmutableSparseNDimArray,sx.MutableSparseNDimArray])
def test_arrays_retain_shape_and_mutability(constructor):
    x=sx.Symbol('x');original=constructor([sx.Rational(1,3),sx.sqrt(x)],(1,2))
    result=sx.nfloat(original,20)
    assert type(result) is type(original)
    assert result.shape==(1,2) and result[0,0]==sx.Float(sx.Rational(1,3),20)
    assert result[0,1]==sx.sqrt(x)


def test_symbolic_matrix_coefficients_float_without_changing_indices():
    x,n=sx.symbols('x n');A=sx.MatrixSymbol('A',n,n)
    result=sx.nfloat(sx.sin(x/2)*A,30)
    assert result==sx.sin(sx.Float(sx.Rational(1,2),30)*x)*A
    assert result.shape==(n,n)
    assert sx.nfloat(A**2,exponent=True)==A**2
    assert sx.nfloat(A[0,1])==A[0,1]


def test_every_nested_power_exponent_is_retained_by_default():
    x,y=sx.symbols('x y')
    value=(sx.sqrt(x)+y)**sx.Rational(2,3)
    assert sx.nfloat(value)==value
    assert sx.nfloat(value,exponent=True).has(sx.Float(sx.Rational(2,3)))


@pytest.mark.parametrize('source',[
    'Lambda(x,0)','Lambda(x,Float(0))','Sum(0,(x,0,3))',
    'Tuple(0,Rational(1,3),sqrt(x))','FiniteSet(0,x/2)',
    'Interval(0,1)','Min(x,Rational(1,3))',
])
def test_zero_and_container_metadata_follow_their_context(source):
    actual=sx.nfloat(expression(sx,source))
    expected=sp.nfloat(expression(sp,source))
    assert actual.to_sympy()==expected


def test_exact_zero_integrals_stay_zero_after_floating_conversion():
    x=sx.Symbol('x')
    result=sx.nfloat(sx.Integral(x,(x,-1,1)))
    assert result.is_zero is True and result.is_Float


def test_undecidable_substitution_and_shared_undefined_arguments_remain_finite():
    x,y=sx.symbols('x y');f,pair=sx.Function('f'),sx.Function('pair')
    value=sx.Subs(sx.Derivative(f(x),x,evaluate=False),x,y)
    assert sx.nfloat(value)==value
    value=sx.Rational(1,3)
    for _ in range(80):value=pair(value,value)
    assert sx.nfloat(value) is value


def test_invalid_precision_and_cyclic_containers_raise():
    with pytest.raises(TypeError):sx.nfloat(1,2.5)
    for digits in [-1,0,10001]:
        assert sx.nfloat(sx.Rational(1,3),digits).to_sympy()==sp.nfloat(sp.Rational(1,3),digits)
    value=[];value.append(value)
    with pytest.raises(ValueError):sx.nfloat(value)


@pytest.mark.parametrize('digits',[-1,0,5,15,40,12001])
def test_float_zero_interchange_retains_numeric_type_and_precision(digits):
    actual=sx.nfloat(0,digits).to_sympy()
    expected=sp.nfloat(0,digits)
    assert type(actual) is type(expected)
    assert actual._prec==expected._prec


def test_nfloat_does_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x=Symbol('x');f=Function('f')
assert nfloat(sqrt(x))==sqrt(x)
assert nfloat(sin(x/2))==sin(Float('0.5')*x)
assert nfloat(f(Rational(1,2)))==f(Rational(1,2))
assert nfloat(x**2,exponent=True)==x**Float(2)
'''
    subprocess.run([sys.executable,'-c',source],check=True)
