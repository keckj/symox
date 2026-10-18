"""Replacement options use native patterns and preserve constructor metadata."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def convert(value):
    if isinstance(value,dict):return {convert(k):convert(v) for k,v in value.items()}
    if isinstance(value,tuple):return tuple(convert(v) for v in value)
    return value.to_sympy() if hasattr(value,'to_sympy') else value


def rules(module):
    x,y=module.symbols('x y');a,b=map(module.Wild,'ab')
    f=module.Function('f')
    coefficient=module.Wild('coefficient',exclude=[x])
    constant=module.Wild('constant',exclude=[x])
    return [
        (module.sin,module.cos),
        (module.sin,lambda arg:module.sin(2*arg)),
        (module.Mul,module.Add),
        (module.sin(a),module.tan(a/2)),
        (module.sin(a),lambda a:module.sin(2*a)),
        (f,lambda arg:arg+1),
        (f(a),lambda a:f(a+1)),
        (coefficient*x+constant,constant-coefficient),
        (lambda item:item.is_Number,lambda item:item**2),
        (lambda item:item.is_Mul,lambda item:2*item),
        (x,y+1),
    ]


@pytest.mark.parametrize('source',[
    'log(sin(x))+tan(sin(x**2))','2*sin(x**3)','x*(x*y+1)',
    '2*x+y','2*x','f(0)+f(1)','sin(sin(x))',
])
@pytest.mark.parametrize('rule',range(11))
def test_replacement_forms_and_maps_match_reference(source,rule):
    expression=sx.sympify(source)
    reference=sp.sympify(source,locals={'f':sp.Function('f')})
    query,value=rules(sx)[rule];Q,V=rules(sp)[rule]
    assert convert(expression.replace(query,value,map=True))==reference.replace(Q,V,map=True)


@pytest.mark.parametrize('exact',[True,False,None])
def test_zero_wildcard_bindings_are_filtered_only_in_exact_mode(exact):
    f=sx.Function('f');F=sp.Function('f');a=sx.Wild('a');A=sp.Wild('a')
    value=f(0)+f(1);reference=F(0)+F(1)
    assert convert(value.replace(f(a),lambda a:f(a+1),exact=exact,map=True)) == \
        reference.replace(F(A),lambda a:F(a+1),exact=exact,map=True)
    assert value.replace(f(a),f(a+1),exact=True)==f(0)+f(2)
    assert value.replace(f(a),f(a+1),exact=False)==f(1)+f(2)


def test_multiple_wildcards_reject_zero_by_default():
    x=sx.Symbol('x');a=sx.Wild('a',exclude=[x]);b=sx.Wild('b',exclude=[x])
    assert (2*x).replace(a*x+b,b-a)==2*x
    assert (2*x).replace(a*x+b,b-a,exact=True)==2*x
    assert (2*x).replace(a*x+b,b-a,exact=False)==-2
    assert (2*x+3).replace(a*x+b,b-a)==1


def test_literal_replacement_targets_complete_nodes():
    x,y=sx.symbols('x y')
    assert (4*x).replace(2*x,y)==4*x
    assert (4*x).subs(2*x,y)==2*y
    held=sx.Add(3*x,4,evaluate=False)
    assert held.replace(3,5)==5*x+4
    assert held.replace(3*x,y)==y+4
    assert held.replace(4,5)==3*x+5


def test_queries_visit_set_children_with_basic_type_flags():
    x=sx.Symbol('x');X=sp.Symbol('x')
    expression=sx.Contains(x,sx.FiniteSet(1,2))
    reference=sp.Contains(X,sp.FiniteSet(1,2))
    assert convert(expression.replace(lambda item:item.is_Number,lambda item:item*item)) == \
        reference.replace(lambda item:item.is_Number,lambda item:item*item)
    assert convert(expression.replace(sx.FiniteSet,lambda *args:sx.FiniteSet(*args,5))) == \
        reference.replace(sp.FiniteSet,lambda *args:sp.FiniteSet(*args,5))


@pytest.mark.parametrize('simultaneous',[True,False])
def test_collapsing_constructor_is_not_rewritten_twice_unless_requested(simultaneous):
    x=sx.Symbol('x');X=sp.Symbol('x')
    expression=sx.sin(sx.sin(x));reference=sp.sin(sp.sin(X))
    actual=expression.replace(sx.sin,lambda arg:x,simultaneous=simultaneous,map=True)
    expected=reference.replace(sp.sin,lambda arg:X,simultaneous=simultaneous,map=True)
    assert convert(actual)==expected
    assert actual[0]==(sx.sin(x) if simultaneous else x)


@pytest.mark.parametrize('source',[
    'Piecewise((x,x>0),(0,True))','Order(x**3,x)',
    'Derivative(f(x),(x,2))','Subs(f(x),x,0)',
    'Integral(x,(x,0,1))','Sum(x,(x,0,3))',
])
def test_replacement_preserves_bound_and_presentation_metadata(source):
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    expression=sx.sympify(source)
    reference=sp.sympify(source,locals={'f':sp.Function('f')})
    assert convert(expression.replace(x,y,map=True))==reference.replace(X,Y,map=True)


def test_undefined_function_head_selection_preserves_assumptions():
    x=sx.Symbol('x');f=sx.Function('f',real=True);g=sx.Function('f',positive=True)
    expression=f(x)+g(x)
    assert expression.replace(f,lambda arg:arg+1)==x+1+g(x)
    fake=sx.Function('sin')
    assert (fake(x)+sx.sin(x)).replace(sx.sin,sx.cos)==fake(x)+sx.cos(x)


def test_predicate_none_return_and_exceptions():
    x=sx.Symbol('x')
    assert sx.sin(x).replace(lambda node:True,lambda node:None)==sx.sin(x)
    def broken(node):raise RuntimeError('callback failed')
    with pytest.raises(RuntimeError,match='callback failed'):x.replace(broken,lambda n:n)
    for query,value in [(sx.sin,1),(lambda e:True,2),(object(),x)]:
        with pytest.raises(TypeError):x.replace(query,value)


def test_compact_shared_graph_is_transformed_without_unfolding():
    x=sx.Symbol('x');expression=x
    for _ in range(75):expression=sx.Mul(expression,expression,evaluate=False)
    seen=[]
    result=expression.replace(lambda node:seen.append(node) or False,lambda node:node)
    assert result is expression
    assert len(seen)==76


def test_matrix_entries_share_the_scalar_replacement_options():
    x=sx.Symbol('x');a=sx.Wild('a')
    matrix=sx.Matrix([[sx.sin(x),2*x],[0,sx.sin(0)]])
    result,mapping=matrix.replace(sx.sin(a),sx.cos(a),map=True,exact=True)
    assert result==sx.Matrix([[sx.cos(x),2*x],[0,0]])
    assert mapping=={sx.sin(x):sx.cos(x)}


def test_scalar_matrix_nodes_preserve_ordered_parents():
    x,y=sx.symbols('x y');A=sx.MatrixSymbol('A',2,2);B=sx.MatrixSymbol('B',2,2)
    expression=sx.Trace(x*A*B)+A[0,1]
    assert expression.replace(x,y)==sx.Trace(y*A*B)+A[0,1]
    assert expression.replace(A,B)==sx.Trace(sx.MatMul(x,B,B))+B[0,1]


def test_replacement_options_do_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x=Symbol('x');a,b=map(Wild,'ab')
assert (2*x).replace(a*x+b,b-a)==2*x
value,mapping=sin(x).replace(sin,cos,map=True)
assert value==cos(x) and mapping=={sin(x):cos(x)}
assert sin(sin(x)).replace(sin,lambda arg:x)==sin(x)
'''
    subprocess.run([sys.executable,'-c',source],check=True)
