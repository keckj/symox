"""Complete hyperbolic solution lattices and principal-branch regressions."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def normalized_lattices(value):
    index=sp.Symbol('_integer_index',integer=True)
    def normalize(image):
        variable=image.lamda.variables[0]
        period=sp.simplify(sp.expand(image.lamda.expr).coeff(variable))
        offset=sp.simplify(image.lamda.expr.subs(variable,0).rewrite(sp.log))
        if period.could_extract_minus_sign():period=-period
        ratio=sp.simplify(sp.re(offset/period))
        if ratio.is_Rational:offset=sp.expand(offset-sp.floor(ratio)*period)
        return sp.ImageSet(sp.Lambda(index,period*index+offset),image.base_set)
    value=value.replace(lambda node:isinstance(node,sp.ImageSet),normalize)
    if isinstance(value,sp.FiniteSet):
        return sp.FiniteSet(*(sp.simplify(item.rewrite(sp.log)) for item in value))
    return value


@pytest.mark.parametrize('head',['sinh','cosh','tanh','coth'])
@pytest.mark.parametrize('target',[-2,-1,0,1,2,sx.I])
@pytest.mark.parametrize('domain',['Reals','Complexes'])
def test_hyperbolic_solutions_match_reference(head,target,domain):
    x=sx.Symbol('x');X=sp.Symbol('x')
    value=getattr(sx,head)(2*x+1)-target
    actual=sx.solveset(value,x,getattr(sx.S,domain)).to_sympy()
    expected=sp.solveset(getattr(sp,head)(2*X+1)-sx.sympify(target).to_sympy(),X,getattr(sp.S,domain))
    assert normalized_lattices(actual)==normalized_lattices(expected)


def test_remote_imaginary_roots_are_not_lost_to_principal_inverse_branches():
    x=sx.Symbol('x')
    odd=sx.solveset(sx.sinh(x),x,sx.S.Complexes)
    even=sx.solveset(sx.cosh(x),x,sx.S.Complexes)
    for n in range(-8,9):
        assert odd.contains(n*sx.I*sx.pi)==sx.true
        assert even.contains((n+sx.Rational(1,2))*sx.I*sx.pi)==sx.true
        assert odd.contains((n+sx.Rational(1,2))*sx.I*sx.pi)==sx.false
    assert sx.solveset(sx.tanh(x)-1,x,sx.S.Complexes)==sx.S.EmptySet
    assert sx.solveset(sx.coth(x)+1,x,sx.S.Complexes)==sx.S.EmptySet


def test_generator_temporary_cannot_capture_a_free_parameter():
    x=sx.Symbol('x');parameter=sx.Symbol('_solveset_t')
    result=sx.solveset(sx.sin(x)+parameter,x,sx.S.Reals)
    assert result.free_symbols=={parameter}
    assert result.subs(parameter,-1).contains(sx.pi/2)==sx.true


def test_unsupported_nonlinear_complex_argument_stays_a_condition_set():
    x=sx.Symbol('x')
    result=sx.solveset(sx.sinh(x*x),x,sx.S.Complexes)
    assert isinstance(result,sx.ConditionSet)


def test_native_hyperbolic_solving_without_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x=Symbol('x')
roots=solveset(sinh(x),x,S.Complexes)
assert roots.contains(7*I*pi)==true
assert solveset(tanh(x)-1,x,S.Complexes)==S.EmptySet
assert solveset(cosh(x),x,S.Reals)==S.EmptySet
'''
    subprocess.run([sys.executable,'-c',source],check=True)
