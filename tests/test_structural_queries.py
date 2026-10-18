"""Find/count selector signatures, public binders and shared DAG multiplicity."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def converted(value):
    if isinstance(value,dict):return {converted(k):converted(v) for k,v in value.items()}
    if isinstance(value,set):return {converted(item) for item in value}
    return value.to_sympy() if hasattr(value,'to_sympy') else value


@pytest.mark.parametrize('head',['Symbol','Integer','Add','Mul','Pow','sin','cos','Function',
    'Expr','Atom','AtomicExpr','Derivative','Integral','Lambda','Order','Tuple','Piecewise','Eq'])
@pytest.mark.parametrize('group',[False,True])
def test_class_and_function_selectors_match_reference(head,group):
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y');f=sx.Function('f');F=sp.Function('f')
    values=[sx.sin(x)+f(x)*sx.cos(x*y),sx.Derivative(f(x),x,2),
        sx.Lambda(x,x+y),sx.Piecewise((x*x,x>0),(y,True)),sx.Order(x*x,x),
        sx.Eq(x,y),sx.S.true,sx.Integral(x*x,x)]
    reference=[sp.sin(X)+F(X)*sp.cos(X*Y),sp.Derivative(F(X),X,2),
        sp.Lambda(X,X+Y),sp.Piecewise((X*X,X>0),(Y,True)),sp.Order(X*X,X),
        sp.Eq(X,Y),sp.S.true,sp.Integral(X*X,X)]
    for native,oracle in zip(values,reference):
        assert converted(native.find(getattr(sx,head),group=group))==oracle.find(getattr(sp,head),group=group)
        assert native.count(getattr(sx,head))==oracle.count(getattr(sp,head))


def test_predicates_literals_wildcards_and_undefined_heads():
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y');f=sx.Function('f');F=sp.Function('f')
    expression=f(x)+sx.sin(x)+sx.cos(y)+x*x
    reference=F(X)+sp.sin(X)+sp.cos(Y)+X*X
    assert converted(expression.find(x,group=True))==reference.find(X,group=True)
    assert converted(expression.find(lambda node:getattr(node,'is_Pow',False)))==reference.find(lambda node:node.is_Pow)
    assert expression.find(f)=={f(x)}
    assert expression.find(F)=={f(x)}
    assert expression.find(sx.Wild('w'))
    collision=sx.Function('sin')
    value=collision(x)+sx.sin(x)
    assert value.find(collision)=={collision(x)}
    assert value.find(sx.sin)=={sx.sin(x)}


def test_matrix_queries_match_undefined_heads_and_preserve_parameter_metadata():
    x=sx.Symbol('x');f=sx.Function('f',real=True)
    matrix=sx.MatrixExpr(sx.Matrix([[f(x),sx.sin(x)],[x,0]]))
    assert matrix.find(f)=={f(x)}
    assert matrix.find(sx.sin)=={sx.sin(x)}
    assert matrix.find(sx.Function('f'))==set()
    assert matrix.count(x)==3


def test_group_queries_count_shared_occurrences_without_unfolding():
    x=sx.Symbol('x');one=sx.S.One
    expression=sx._core._restore_arithmetic('Mul',sx.Integer(2),[(x,one)],False)
    for _ in range(80):
        expression=sx._core._restore_arithmetic('Mul',one,[(expression,one),(expression,one)],False)
    assert expression.find(x,group=True)=={x:2**80}
    assert expression.count(x)==2**80
    assert expression.find(sx.Symbol)=={x}


def test_native_node_identity_preserves_alpha_equal_binder_presentations():
    x,y=sx.symbols('x y');f=sx.Function('f')
    left=sx.Subs(f(x),x,0);right=sx.Subs(f(y),y,0)
    assert left==right and left._dag_id!=right._dag_id
    expression=sx.Add(left,right,evaluate=False)
    assert expression.find(x)=={x}
    assert expression.find(y)=={y}
    assert expression.find(sx.Symbol,group=True)=={x:2,y:2}
    assert expression.has(y.func)


def test_foreign_constructor_selectors_do_not_change_native_outputs():
    x=sx.Symbol('x');expression=sx.sin(x)+x*x
    assert expression.find(sp.sin)=={sx.sin(x)}
    assert expression.find(sp.Symbol)=={x}


def test_without_sympy_import():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x=Symbol('x');f=Function('f')
expression=sin(x)+f(x)
assert expression.find(sin)=={sin(x)}
assert expression.find(Symbol,group=True)=={x:2}
assert expression.find(f)=={f(x)}
'''
    subprocess.run([sys.executable,'-c',source],check=True)
