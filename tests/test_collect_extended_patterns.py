"""Discovered collection generators and Wild patterns against the reference."""
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp
from test_collect_compat import CASES,converted,native,oracle


@pytest.mark.parametrize('source,patterns',CASES+[
    ('x*exp(x)+sin(x)*y+2*sin(x)+3*x',['x']),
    ('a*x*y+x*y+b*x+x',['x','y']),
    ('a*x**2+b*x**2+c*x**3+d*x**3',['x']),
    ('a*sqrt(x)+b*sqrt(x)+c*x+d*x',['x']),
    ('a*sin(x)*x+b*sin(x)*x+c*sin(x)',['x']),
    ('a*x**y+b*x**y+c*x**(2*y)',['x']),
    ('a*Integral(x,(x,0,1))+b*Integral(x,(x,0,1))',['x']),
])
@pytest.mark.parametrize('evaluate',[False,True])
def test_exact_none_discovers_generators_and_retains_partitions(source,patterns,evaluate):
    expression=native(source)
    actual=sx.collect(expression,[native(p) for p in patterns],exact=None,evaluate=evaluate)
    expected=sp.collect(oracle(source),[oracle(p) for p in patterns],exact=None,evaluate=evaluate)
    assert converted(actual)==expected
    reconstructed=sum((k*v for k,v in actual.items()),sx.Integer(0)) if isinstance(actual,dict) else actual
    assert sx.expand(reconstructed-expression)==0


WILD_CASES=[
    ('a*x**y-b*x**y','w**y'),
    ('a*x**(3*y)+b*x**(3*y)+c*x**y','w**y'),
    ('a*x**y+b*x**y+c*z**y','w**v'),
    ('a*sin(x)+b*sin(x)+c*cos(x)','sin(w)'),
    ('a*f(x)+b*f(x)+c*f(y)','f(w)'),
    ('a*x**2+b*x**2+c*y','w'),
    ('a*x+b*x+c*y','w*x'),
    ('a*x**2+b*x**2+c*x**3','x**w'),
    ('a*x**y+b*x**y+c*z**y','w**y'),
]


def wild_expression(module,source,pattern,exclude):
    w=module.Wild('w',exclude=[module.Symbol('z')] if exclude else [])
    v=module.Wild('v')
    return module.sympify(source),module.sympify(pattern,locals={'w':w,'v':v})


@pytest.mark.parametrize('source,pattern',WILD_CASES)
@pytest.mark.parametrize('exact',[False,True,None])
@pytest.mark.parametrize('evaluate',[False,True])
@pytest.mark.parametrize('exclude',[False,True])
def test_wild_patterns_keep_actual_factor_keys(source,pattern,exact,evaluate,exclude):
    expression,p=wild_expression(sx,source,pattern,exclude)
    reference,P=wild_expression(sp,source,pattern,exclude)
    actual=sx.collect(expression,p,exact=exact,evaluate=evaluate)
    expected=sp.collect(reference,P,exact=exact,evaluate=evaluate)
    assert converted(actual)==expected
    reconstructed=sum((k*v for k,v in actual.items()),sx.Integer(0)) if isinstance(actual,dict) else actual
    assert sx.expand(reconstructed-expression)==0


@pytest.mark.parametrize('evaluate',[False,True])
def test_discovery_keeps_callback_coefficient_semantics(evaluate):
    source='a*x*y+x*y+b*x+x+c*sin(x)+d*sin(x)'
    actual=sx.collect(native(source),native('x'),func=lambda c:c+1,exact=None,evaluate=evaluate)
    expected=sp.collect(oracle(source),oracle('x'),func=lambda c:c+1,exact=None,evaluate=evaluate)
    assert converted(actual)==expected


def test_extended_collection_never_imports_reference():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('reference import')
sys.meta_path.insert(0,Block())
from symox import *
x,y,a,b=symbols('x y a b');w=Wild('w')
assert collect(a*x**y-b*x**y,w**y)==x**y*(a-b)
assert collect(a*x*y+x*y+b*x+x,[x,y],exact=None)==x*y*(a+1)+x*(b+1)
'''
    subprocess.run([sys.executable,'-c',source],check=True)
