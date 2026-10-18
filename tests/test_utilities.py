"""Container oracles, traversal behavior, graph invariants and memoization."""
import itertools
import random
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from symox.utilities import (sift,flatten,unflatten,reshape,group,subsets,variations,
    cartes,prefixes,postfixes,rotations,has_dups,has_variety,take,prod,dict_merge,
    topological_sort,ordered,preorder_traversal,postorder_traversal,bottom_up,use,
    cacheit,memoize_property)


def normalized(value):
    if isinstance(value,(list,tuple)):
        return type(value)(normalized(item) for item in value)
    return value.to_sympy() if hasattr(value,'to_sympy') else value


@pytest.mark.parametrize('name',['prefixes','postfixes','rotations','group','has_dups','has_variety'])
@pytest.mark.parametrize('seq',[[],[1],[1,1,2,3,3],(1,2,3),'aabcc'])
def test_sequence_helpers_match_reference(name,seq):
    actual=globals()[name](seq);expected=getattr(sp,name)(seq)
    if name in ('prefixes','postfixes','rotations'):
        actual,expected=list(actual),list(expected)
    assert actual==expected


@pytest.mark.parametrize('length',range(5))
@pytest.mark.parametrize('repetition',[False,True])
def test_combinatorial_iterators(length,repetition):
    for size in range(5):
        assert list(subsets(range(length),size,repetition))==list(sp.subsets(range(length),size,repetition))
        assert list(variations(range(length),size,repetition))==list(sp.variations(range(length),size,repetition))
    assert list(subsets(range(length),repetition=repetition))==list(sp.subsets(range(length),repetition=repetition))
    assert list(cartes(range(length),'ab'))==list(itertools.product(range(length),'ab'))


def test_flatten_depth_symbolic_class_and_array_inputs():
    nested=[1,[2,(3,[4])],5]
    for levels in [None,0,1,2,3]:
        assert flatten(nested,levels)==sp.flatten(nested,levels)
    x,y,z=sx.symbols('x y z');X,Y,Z=sp.symbols('x y z')
    assert normalized(flatten([x+y,z],cls=sx.Add))==sp.flatten([X+Y,Z],cls=sp.Add)
    assert flatten([sx.Array([[1,2],[3,4]])])==[1,2,3,4]
    iterator=iter([1,2])
    assert flatten([iterator])==[iterator]
    assert flatten(nested,0) is nested
    cyclic=[];cyclic.append(cyclic)
    with pytest.raises(ValueError):flatten(cyclic)
    with pytest.raises(ValueError):flatten(nested,-1)


@pytest.mark.parametrize('template',[[4],(4,),(2,2),(2,[2]),((2,),[2]),(1,[2],1),
    ([[1],1,(2,)],),([1],1,(2,))])
@pytest.mark.parametrize('kind',[list,tuple])
def test_nested_reshape_templates_match_reference(template,kind):
    values=kind(range(8))
    assert reshape(values,template)==sp.reshape(values,template)
    assert sorted(flatten(reshape(values,template)))==list(range(8))


def test_grouping_sifting_and_streaming_helpers():
    assert group([1,1,2,3,3],False)==[(1,2),(2,1),(3,2)]
    assert unflatten(iter(range(6)),3)==[(0,1,2),(3,4,5)]
    for size in [0,-1,4,1.5]:
        with pytest.raises(ValueError):unflatten(range(6),size)
    assert sift(range(6),lambda value:value%2)==sp.sift(range(6),lambda value:value%2)
    groups=sift([1],lambda value:value)
    assert groups[99]==[]
    assert sift(range(4),lambda value:value%2,binary=True)==([1,3],[0,2])
    with pytest.raises(ValueError):sift([0],lambda value:2,binary=True)
    assert take(itertools.count(),4)==[0,1,2,3]
    assert prod(range(1,6),2)==240
    assert type(prod([2,3])) is int
    assert dict_merge({'a':1},{'b':2,'a':3})=={'a':3,'b':2}
    assert has_dups([[1],[2],[1]]) and not has_dups(iter([[1],[2]]))


@pytest.mark.parametrize('value',['x',b'ab',[1,'x'],(1,[2,3]),{3,1,2},
    {'x':(1,2),'y':[3]},None,True,3.5])
def test_generic_sort_keys_match_reference(value):
    assert normalized(sx.default_sort_key(value))==sp.default_sort_key(value)


@pytest.mark.parametrize('kind',['ImmutableDenseNDimArray','MutableDenseNDimArray',
    'ImmutableSparseNDimArray','MutableSparseNDimArray'])
def test_array_sort_keys_and_sparse_resource_bounds(kind):
    native=getattr(sx,kind)([[1,0],[2,3]])
    reference=getattr(sp,kind)([[1,0],[2,3]])
    assert normalized(sx.default_sort_key(native))==sp.default_sort_key(reference)
    large=sx.ImmutableSparseNDimArray({5:2},(10**12,))
    assert sx.default_sort_key(large)[0]==(5,0,'ImmutableSparseNDimArray')
    with pytest.raises(ValueError):sx.default_sort_key(large.as_mutable())


def test_ordered_default_and_lazy_tie_breaking():
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    native=[x+y,x**2,x,y,sx.sin(x),[x,y],[1,x],(x,)]
    reference=[X+Y,X**2,X,Y,sp.sin(X),[X,Y],[1,X],(X,)]
    assert normalized(list(ordered(native)))==list(sp.ordered(reference))
    calls=[]
    def second(value):calls.append(value);return value
    iterator=ordered([3,1,2],keys=[lambda value:value,second],default=False)
    assert list(iterator)==[1,2,3] and calls==[]
    assert list(ordered([2,1],keys=lambda value:0,default=False))==[2,1]
    with pytest.raises(ValueError):list(ordered([2,1],keys=lambda value:0,default=False,warn=True))
    with pytest.raises(ValueError):list(ordered([1],default=False))


@pytest.mark.parametrize('traversal',[preorder_traversal,postorder_traversal])
def test_traversal_matches_scalar_tree_and_python_containers(traversal):
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    expression=(x+y)**2+sx.sin(x*y)
    reference=(X+Y)**2+sp.sin(X*Y)
    oracle=getattr(sp,traversal.__name__)
    assert normalized(list(traversal(expression)))==list(oracle(reference))
    assert list(traversal([1,[2,3]]))==list(oracle([1,[2,3]]))
    assert normalized(list(traversal(expression,keys=True)))==list(oracle(reference,keys=True))


def test_preorder_skip_and_deep_tree_iteration():
    x,y=sx.symbols('x y');expression=(x+y)**2
    traversal=preorder_traversal(expression)
    assert next(traversal)==expression
    assert next(traversal)==x+y
    traversal.skip()
    assert list(traversal)==[sx.Integer(2)]
    nested=1
    for _ in range(2000):nested=[nested]
    assert sum(1 for _ in preorder_traversal(nested))==2001
    assert sum(1 for _ in postorder_traversal(nested))==2001


def test_tree_transformations_respect_levels_and_atoms():
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    expr=(x+y)**2*x+1;reference=(X+Y)**2*X+1
    for level in range(5):
        assert normalized(use(expr,sx.expand,level))==sp.use(reference,sp.expand,level)
    assert normalized(bottom_up(expr,sx.expand))==sp.bottom_up(reference,sp.expand)
    def replace(value):return y if value==x else value
    assert bottom_up(sx.sin(x),replace,atoms=True)==sx.sin(y)
    assert bottom_up(sx.sin(x),replace)==sx.sin(x)
    assert bottom_up(2,lambda value:value+1,nonbasic=True)==3


def test_topological_order_oracle_and_graph_invariants():
    vertices=[2,3,5,7,8,9,10,11]
    edges=[(7,11),(7,8),(5,11),(3,8),(3,10),(11,2),(11,9),(11,10),(8,9)]
    for key in [None,lambda value:-value]:
        assert topological_sort((vertices,edges),key)==sp.topological_sort((vertices,edges),key)
    rng=random.Random(203)
    for count in range(1,30):
        edges=[(i,j) for i in range(count) for j in range(i+1,count) if rng.random()<.2]
        ordering=topological_sort((list(range(count)),edges))
        positions={value:i for i,value in enumerate(ordering)}
        assert set(ordering)==set(range(count))
        assert all(positions[a]<positions[b] for a,b in edges)
    with pytest.raises(ValueError):topological_sort(([1,2],[(1,2),(2,1)]))
    with pytest.raises(ValueError):topological_sort(([1],[(1,2)]))


def test_memoization_types_unhashable_inputs_exceptions_and_invalidation():
    calls=[]
    @cacheit
    def value(item):calls.append(item);return item
    assert value(1)==value(1)==1 and len(calls)==1
    assert value(1.0)==1 and len(calls)==2
    assert value([1])==value([1])==[1] and len(calls)==4
    assert value.cache_info().hits==1
    value.cache_clear()
    value(1);assert len(calls)==5
    errors=[]
    @cacheit
    def fails(item):errors.append(item);raise TypeError('unhashable type: from inside function')
    with pytest.raises(TypeError):fails(1)
    assert errors==[1]
    class WithProperty:
        __slots__=('_answer','calls')
        def __init__(self):self.calls=0
        @memoize_property
        def answer(self):self.calls+=1;return 42
    instance=WithProperty()
    assert instance.answer==instance.answer==42 and instance.calls==1


def test_no_sympy_import():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import symbols,sin
from symox import ordered,postorder_traversal,flatten,topological_sort,cacheit
x,y=symbols('x y')
assert list(ordered([y,x]))==[x,y]
assert list(postorder_traversal(sin(x)))==[x,sin(x)]
assert flatten([1,[2,[3]]])==[1,2,3]
assert topological_sort(([1,2],[(1,2)]))==[1,2]
'''
    subprocess.run([sys.executable,'-c',source],check=True)
