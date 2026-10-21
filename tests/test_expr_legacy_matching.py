"""Legacy scalar matching can solve a linear wildcard after native matching."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('source,pattern',[
    ('2/x','a*x'),('2/(x+1)','a*(x+1)'),('sin(x)','a*x'),
    ('2','a*x'),('1/x','a*x'),('f(0)','a*x'),
    ('x-2','a-x'),('x+1','a-x'),('x+y','a*x'),
    ('2*x+3','a*x+3'),('x*x+1','a*x*x+1'),
])
def test_linear_legacy_matching_matches_reference(source,pattern):
    x=sx.Symbol('x');a=sx.Wild('a');A=sp.Wild('a')
    expression=sx.sympify(source);reference=sp.sympify(source,locals={'f':sp.Function('f')})
    query=sx.sympify(pattern).subs(sx.Symbol('a'),a)
    Q=sp.sympify(pattern,locals={'a':A})
    actual=expression.match(query,old=True)
    expected=reference.match(Q,old=True)
    assert (actual is None)==(expected is None)
    if actual is not None:
        assert {key.to_sympy():value.to_sympy() for key,value in actual.items()}==expected
        assert sx.expand(query.xreplace(actual)-expression)==0
    assert query.matches(expression,old=True)==actual


def test_default_flag_and_exclusions_preserve_native_matching():
    x,y=sx.symbols('x y');a=sx.Wild('a')
    assert (2/x).match(a*x) is None
    assert (2/x).match(a*x,old=False) is None
    assert (2/x).match(a*x,old=True)=={a:2/x**2}
    forbidden=sx.Wild('a',exclude=[x])
    assert (2/x).match(forbidden*x,old=True) is None
    assert y.match(a*x,old=True) is None
    assert x.match(a*x,old=True)=={a:sx.S.One}


def test_preexisting_bindings_are_copied_in_function_argument_matching():
    x=sx.Symbol('x');a,b=map(sx.Wild,'ab');f=sx.Function('f')
    bindings={a:x}
    query=f(a,b)
    assert query.matches(f(x,2),bindings,old=True)=={a:x,b:sx.Integer(2)}
    assert bindings=={a:x}
    assert query.matches(f(x+1,2),bindings,old=True) is None


def test_legacy_fallback_matches_the_reference_on_nonlinear_patterns():
    a=sx.Wild('a')
    # SymPy's Mul.matches turns a*a into a**2 and Pow.matches takes the
    # principal root.
    assert sx.Integer(2).match(a*a,old=True)=={a:sx.sqrt(2)}
    assert sx.Integer(2).match(sx.sin(a),old=True) is None


def test_matrix_matching_can_use_legacy_scalar_coefficient_bindings():
    x=sx.Symbol('x');a=sx.Wild('a');A=sx.MatrixSymbol('A',2,2);B=sx.MatrixSymbol('B',2,2)
    expression=(2/x)*A*B;pattern=(a*x)*A*B
    assert expression.match(pattern) is None
    assert expression.match(pattern,old=True)=={a:2/x**2}
    assert pattern.matches(expression,old=True)=={a:2/x**2}
    assert expression.match((a*x)*B*A,old=True) is None


def test_legacy_matching_does_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x=Symbol('x');a=Wild('a')
assert (2/x).match(a*x,old=True)=={a:2/x**2}
assert (a*x).matches(2/x,old=True)=={a:2/x**2}
assert (2/x).match(a*x,old=False) is None
'''
    subprocess.run([sys.executable,'-c',source],check=True)
