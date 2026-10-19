"""Callable queries pass the whole matched expression to their replacement."""
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp
from symox.codegen.rewriting import ReplaceOptim,exp2_opt,log1p_opt
from sympy.codegen.rewriting import ReplaceOptim as Reference


@pytest.mark.parametrize('source',['x**2','sin(x)**3','x**2+y**3','sin(x**2)','(x+y)**3'])
def test_callable_replacement_receives_one_expression(source):
    value=sx.sympify(source);expected=sp.sympify(source);seen=[]
    def replace(expression):
        seen.append(expression)
        return expression.base+expression.exp
    actual=ReplaceOptim(lambda node:node.is_Pow,replace)(expr=value)
    reference=Reference(lambda node:node.is_Pow,lambda node:node.base+node.exp)(expr=expected)
    assert actual.to_sympy()==reference
    assert seen and all(node.is_Pow for node in seen)


def test_function_class_query_still_passes_function_arguments():
    x=sx.Symbol('x');X=sp.Symbol('x')
    actual=ReplaceOptim(sx.sin,sx.cos)(sx.sin(x)+sx.sin(x*x))
    expected=Reference(sp.sin,sp.cos)(sp.sin(X)+sp.sin(X*X))
    assert actual.to_sympy()==expected


def test_atom_queries_are_not_skipped():
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    for expression,reference in [(x,X),(x+1,X+1)]:
        actual=ReplaceOptim(lambda node:node==x,lambda node:y)(expression)
        expected=Reference(lambda node:node==X,lambda node:Y)(reference)
        assert actual.to_sympy()==expected


def test_builtin_rewrites_keep_their_correct_arity():
    x=sx.Symbol('x')
    assert str(exp2_opt(2**x))=='exp2(x)'
    assert str(log1p_opt(sx.log(1+x)))=='log1p(x)'


def test_rewrite_callbacks_do_not_import_reference():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('reference import')
sys.meta_path.insert(0,Block())
from symox import *
from symox.codegen.rewriting import ReplaceOptim,exp2_opt
x=Symbol('x')
assert ReplaceOptim(lambda n:n.is_Pow,lambda n:n.base+n.exp)(x*x)==x+2
assert str(exp2_opt(2**x))=='exp2(x)'
'''
    subprocess.run([sys.executable,'-c',source],check=True)
