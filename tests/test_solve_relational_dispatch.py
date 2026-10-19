"""Relational solve results retain every supplied inequality and equation."""
import inspect
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize('source',[
    'x<3','x**2>4','Ne(x,1)',
    '[x<3,x**2>4]','[x>3,x<1]',
    '[x>0,x**2-1]','[x>=0,Eq(x**2,4)]',
    '[1/(x-1)>0,x<4]','[Abs(x-1)<2,x>0]',
])
@pytest.mark.parametrize('flags',[{}, {'dict':True}, {'set':True}])
def test_relational_solution_sets_match_reference(source,flags):
    def build(module):
        x=module.Symbol('x')
        return eval(source,dict(vars(module),x=x)),x
    expression,x=build(sx);reference,X=build(sp)
    actual=sx.solve(expression,x,**flags)
    expected=sp.solve(reference,X,**flags)
    assert actual.to_sympy().as_set()==expected.as_set()


@pytest.mark.parametrize('real',[False,True])
def test_parameter_equation_is_kept_together_with_the_inequality(real):
    def build(module):
        x,y=module.symbols('x y',**({'real':True} if real else {}))
        return [x+y-3,x>3],x
    expression,x=build(sx);reference,X=build(sp)
    actual=sx.solve(expression,x).to_sympy()
    expected=sp.solve(reference,X)
    assert sp.simplify(sp.Equivalent(actual,expected)) is sp.true


@pytest.mark.parametrize('container',['scalar','list','tuple','native_tuple','iterator','string'])
def test_inferred_variable_and_input_sequences_keep_relational_output(container):
    x=sx.Symbol('x');expression=x<3
    inputs={'scalar':expression,'list':[expression],'tuple':(expression,),
            'native_tuple':sx.Tuple(expression),'iterator':iter([expression]),'string':'x<3'}
    actual=sx.solve(inputs[container])
    X=sp.Symbol('x')
    assert actual.to_sympy().as_set()==sp.solve(X<3).as_set()


def test_public_alias_signature_and_algebraic_outputs_are_preserved():
    from symox.solvers.solvers import solve
    assert solve is sx.solve
    assert str(inspect.signature(solve))==str(inspect.signature(sp.solve))
    x=sx.Symbol('x')
    assert solve(x-2,x)==[sx.Integer(2)]
    assert solve([x-2],[x])=={x:sx.Integer(2)}


def test_relational_dispatch_does_not_import_reference():
    script='''
import sys
class NoReference:
    def find_spec(self,fullname,path=None,target=None):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise ImportError('reference forbidden')
sys.meta_path.insert(0,NoReference())
import symox as s
from symox.solvers.solvers import solve
x,y=s.symbols('x y')
assert solve is s.solve
assert solve(x<3,x)==s.reduce_inequalities([x<3],x)
assert solve([x+y-3,x>3],x)==s.reduce_inequalities([x+y-3,x>3],x)
'''
    subprocess.run([sys.executable,'-c',script],check=True,timeout=20)
