"""Scalar substitution semantics are shared by all concrete matrix entries."""
from collections import UserDict
from types import MappingProxyType
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def pairs(x,y,container):
    substitutions=[(x,y),(y,2)]
    if container=='mapping':return {x:y,y:2}
    if container=='userdict':return UserDict(substitutions)
    if container=='proxy':return MappingProxyType(dict(substitutions))
    if container=='list':return [[old,new] for old,new in substitutions]
    if container=='tuple':return tuple(substitutions)
    if container=='generator':return iter(substitutions)
    if container=='nested_generator':return (iter(item) for item in substitutions)
    raise AssertionError(container)


@pytest.mark.parametrize('kind',['Matrix','ImmutableMatrix','SparseMatrix','ImmutableSparseMatrix'])
@pytest.mark.parametrize('container',['mapping','userdict','proxy','list','tuple','generator','nested_generator'])
@pytest.mark.parametrize('simultaneous',[False,True])
def test_shared_substitutions_match_reference(kind,container,simultaneous):
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    matrix=getattr(sx,kind)([[x+y,x/y],[x*x+x*y,sx.sin(x)]])
    reference=getattr(sp,kind)([[X+Y,X/Y],[X*X+X*Y,sp.sin(X)]])
    original=matrix.copy()
    actual=matrix.subs(pairs(x,y,container),simultaneous=simultaneous)
    # SymPy matrices reuse already exhausted inner iterators. Compare the
    # same rules materialized once; the native adapter accepts both levels.
    expected=reference.subs(pairs(X,Y,'tuple' if container=='nested_generator' else container),simultaneous=simultaneous)
    assert actual.to_sympy()==expected and matrix==original
    assert type(actual) is type(matrix)
    assert isinstance(actual,sx.ImmutableMatrix)==isinstance(matrix,sx.ImmutableMatrix)


@pytest.mark.parametrize('kind',['Matrix','ImmutableMatrix','SparseMatrix','ImmutableSparseMatrix'])
@pytest.mark.parametrize('simultaneous',[False,True])
def test_symbol_swaps_and_simultaneous_zero_denominators(kind,simultaneous):
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    matrix=getattr(sx,kind)([[x,y],[x/y,0]])
    reference=getattr(sp,kind)([[X,Y],[X/Y,0]])
    for rules,foreign in [({x:y,y:x},{X:Y,Y:X}),({x:0,y:0},{X:0,Y:0})]:
        actual=matrix.subs(rules,simultaneous=simultaneous)
        expected=reference.subs(foreign,simultaneous=simultaneous)
        assert actual.to_sympy()==expected


@pytest.mark.parametrize('kind',['Matrix','ImmutableMatrix','SparseMatrix','ImmutableSparseMatrix'])
def test_one_shot_substitutions_also_work_without_options(kind):
    x,y=sx.symbols('x y')
    matrix=getattr(sx,kind)([[x,x],[x,x]])
    assert matrix.subs(iter([(x,y),(y,2)]))==2*sx.ones(2)
    assert matrix.subs(iter([(0,1)]))==matrix
    assert getattr(sx,kind)([[0,x]]).subs(iter([(0,1)]))==sx.Matrix([[1,x]])


@pytest.mark.parametrize('kind',['Matrix','ImmutableMatrix','SparseMatrix','ImmutableSparseMatrix'])
@pytest.mark.parametrize('shape',[(0,0),(0,3),(2,0)])
def test_empty_matrix_substitution_has_no_entries_to_visit(kind,shape):
    matrix=getattr(sx,kind)(*shape,[])
    assert matrix.subs(simultaneous=True)==matrix


@pytest.mark.parametrize('container',['userdict','proxy','list','tuple','generator','nested_generator'])
def test_scalar_substitution_accepts_mapping_and_pair_protocols(container):
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    for simultaneous in [False,True]:
        assert (x+y).subs(pairs(x,y,container),simultaneous=simultaneous).to_sympy()==(X+Y).subs(pairs(X,Y,container),simultaneous=simultaneous)


@pytest.mark.parametrize('rules',[[],[1],[1,2,3]])
def test_malformed_pair_lengths_are_rejected(rules):
    x=sx.Symbol('x')
    for expression in [x,sx.Matrix([[x]]),sx.SparseMatrix([[x]])]:
        with pytest.raises(ValueError):expression.subs([rules])


def test_mapping_access_and_iterator_errors_keep_the_original_exception():
    class Marker(Exception):pass
    failure=Marker('mapping failed')
    class BadMapping:
        @property
        def items(self):raise failure
        def __iter__(self):return iter([])
    x=sx.Symbol('x')
    with pytest.raises(Marker) as caught:x.subs(BadMapping())
    assert caught.value is failure
    def bad():
        yield x,1
        raise failure
    with pytest.raises(Marker) as caught:sx.Matrix([[x]]).subs(bad())
    assert caught.value is failure


@pytest.mark.parametrize('kind',['Matrix','ImmutableMatrix','SparseMatrix','ImmutableSparseMatrix'])
@pytest.mark.parametrize('simultaneous',[False,True])
def test_held_derivative_substitution_retains_scalar_scope_rules(kind,simultaneous):
    x,y=sx.symbols('x y');f=sx.Function('f')
    derivative=sx.Derivative(f(x),x,evaluate=False)
    matrix=getattr(sx,kind)([[derivative,0]])
    result=matrix.subs(x,y*y,simultaneous=simultaneous)
    assert result[0,0]==derivative.subs(x,y*y,simultaneous=simultaneous)
    assert result[0,1]==0


def test_mapping_priority_checks_do_not_expand_shared_or_high_order_patterns():
    x,y,z=sx.symbols('x y z');f=sx.Function('f')
    derivative=sx.Derivative(f(x),(x,10**9),evaluate=False)
    assert derivative.subs({derivative:y,x:z})==y
    expression=x
    for _ in range(80):expression=f(expression,expression)
    assert expression.subs({expression:y,x:z})==y


def test_matrix_substitution_does_not_import_sympy():
    source='''
import sys
from collections import UserDict
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x,y=symbols('x y')
for cls in [Matrix,ImmutableMatrix,SparseMatrix,ImmutableSparseMatrix]:
    matrix=cls([[x,y],[x/y,0]])
    assert matrix.subs(UserDict({x:y,y:x}),simultaneous=True)==cls([[y,x],[y/x,0]])
    assert matrix.subs(iter([(x,y),(y,2)]))==cls([[2,2],[1,0]])
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=30)
