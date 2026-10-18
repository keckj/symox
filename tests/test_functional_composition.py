"""Functional chains must reconstruct their original native expression."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from symox.solvers.decompogen import compogen
from sympy.solvers.decompogen import compogen as reference_compogen


@pytest.mark.parametrize('source',[
    'sin(cos(x))','sin(x)**2+sin(x)+1','sqrt(6*x*x-5)',
    'sin(sqrt(cos(x*x+1)))','x**4+2*x**3-x-1',
    'x**4+2*x*x+1','(x*x+1)**3','exp(x*x+1)',
    'log(sin(x))','Abs(x*x-1)','1/(x*x+1)',
    '(sin(x)+1)/(sin(x)-1)','sin(x)+cos(x)',
    'x**x','x**(x+1)','x','2','a*x*x+1',
    'Min(sin(x)**2+1,sin(x)**3)',
    'Max(cos(x)**2,2*cos(x)+1,3)',
    'Min(sin(x),cos(x))','f(sin(x))','f(x,x+1)',
])
def test_general_decompositions_match_reference_and_reconstruct(source):
    x=sx.Symbol('x');X=sp.Symbol('x')
    expression=sx.sympify(source);foreign=sp.sympify(source)
    chain=sx.decompogen(expression,x)
    expected=sp.decompogen(foreign,X)
    assert [part.to_sympy() for part in chain]==expected
    assert sx.simplify(compogen(chain,x)-expression)==0
    assert sp.simplify(reference_compogen(expected,X)-foreign)==0


def test_variable_exponents_keep_the_base_and_reconstruct():
    x=sx.Symbol('x')
    expression=2**(x*x+1)
    assert sx.decompogen(expression,x)==[2**x,x*x+1]
    assert compogen(sx.decompogen(expression,x),x)==expression
    expression=(x+1)**(x*x+1)
    assert sx.decompogen(expression,x)==[expression]


def test_parameters_and_multiple_argument_functions_do_not_get_captured():
    x,a=sx.symbols('x a');f=sx.Function('f')
    expression=f(sx.sin(x),a)
    assert sx.decompogen(expression,x)==[f(x,a),sx.sin(x)]
    assert compogen(sx.decompogen(expression,x),x)==expression
    assert sx.decompogen(a,x)==[a]
    assert compogen((part for part in [sx.sin(x),x+a]),x)==sx.sin(x+a)


@pytest.mark.parametrize('value',[sx.true,sx.Eq(sx.Symbol('x'),0),sx.FiniteSet(1),sx.Matrix([[1]])])
def test_only_arithmetic_scalar_functions_are_composed(value):
    x=sx.Symbol('x')
    with pytest.raises(TypeError):sx.decompogen(value,x)
    with pytest.raises(TypeError):compogen([value],x)


def test_invalid_variables_and_empty_or_excessive_chains_raise():
    x=sx.Symbol('x')
    with pytest.raises(ValueError):sx.decompogen(x,x*x)
    with pytest.raises(ValueError):compogen([],x)
    with pytest.raises(NotImplementedError):compogen([x]*257,x)


def test_composition_uses_no_sympy_runtime():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
from symox.solvers.decompogen import compogen
x=Symbol('x')
expression=sin(sqrt(cos(x*x+1)))
chain=decompogen(expression,x)
assert chain==[sin(x),sqrt(x),cos(x),x*x+1]
assert compogen(chain,x)==expression
assert decompogen(x**x,x)==[x**x]
'''
    subprocess.run([sys.executable,'-c',source],check=True)
