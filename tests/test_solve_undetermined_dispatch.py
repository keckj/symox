"""Scalar solve attempts coefficient identities before algebraic families."""
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp


def foreign(value):
    if isinstance(value,dict):return {foreign(k):foreign(v) for k,v in value.items()}
    if isinstance(value,(list,tuple,set)):return type(value)(foreign(v) for v in value)
    return value.to_sympy() if hasattr(value,'to_sympy') else value


@pytest.mark.parametrize('source', [
    'a*x+b-2*x-3', '(a*x+b-2*x-3)/(x+1)',
    'a*x**100+b-2*x**100-3', 'a*x+b*y-2*x-3*y', '(a-b)*x',
])
@pytest.mark.parametrize('ordering', ['positional','list','reversed','set'])
@pytest.mark.parametrize('flags', [{},{'dict':True},{'set':True}])
def test_unique_coefficient_solution_and_output_conventions(source, ordering, flags):
    x,y,a,b=sx.symbols('x y a b')
    expression=eval(source,dict(x=x,y=y,a=a,b=b))
    variables={'positional':(a,b),'list':([a,b],),'reversed':(b,a),'set':({a,b},)}[ordering]
    assert foreign(sx.solve(expression,*variables,**flags))==sp.solve(
        expression.to_sympy(),*foreign(variables),**flags)


def test_a_list_requests_an_algebraic_system_instead_of_an_identity():
    x,a,b=sx.symbols('x a b');expression=a*x+b-2*x-3
    assert sx.solve(expression,a,b)=={a:2,b:3}
    ordinary=sx.solve([expression],a,b)
    reference=sp.solve([expression.to_sympy()],a.to_sympy(),b.to_sympy())
    assert foreign(ordinary)==reference


def test_an_impossible_coefficient_identity_falls_back_to_algebraic_solving():
    x,a,b=sx.symbols('x a b');expression=a*x+b-x*x
    result=sx.solve(expression,a,b)
    assert isinstance(result,list) and len(result)==1
    row=result[0]
    assert len(row)==2
    assert sx.cancel(expression.subs(dict(zip((a,b),row))))==0


def test_sparse_high_degree_identity_dispatch_never_imports_reference():
    script='''
import sys
class NoReference:
    def find_spec(self, fullname, path=None, target=None):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise ImportError('reference forbidden')
sys.meta_path.insert(0,NoReference())
import symox as sx
x,a,b=sx.symbols('x a b')
assert sx.solve(a*x**1000000+b-2*x**1000000-3,a,b,cubics=False)=={a:2,b:3}
assert len(sx.roots(x**3-x+1,x))==3
assert sx.solve([x*x-2,a-3],x)==[(-sx.sqrt(2),),(sx.sqrt(2),)]
'''
    subprocess.run([sys.executable,'-c',script],check=True,timeout=30)
