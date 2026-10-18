"""Comprehensions expand native expressions and preserve their bound indices."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from sympy.tensor.array import ArrayComprehension as Reference


@pytest.mark.parametrize('formula',['10*i+j','i*j','sin(i)+j','i**2+a*j'])
@pytest.mark.parametrize('bounds',[((1,3),(1,2)),((-2,1),(2,4)),((0,0),(0,0))])
def test_expansion_and_metadata_match_reference(formula,bounds):
    i,j,a=sx.symbols('i j a');I,J,A=sp.symbols('i j a')
    f=eval(formula,vars(sx)|{'i':i,'j':j,'a':a})
    g=eval(formula,vars(sp)|{'i':I,'j':J,'a':A})
    actual=sx.ArrayComprehension(f,(i,*bounds[0]),(j,*bounds[1]))
    expected=Reference(g,(I,*bounds[0]),(J,*bounds[1]))
    assert actual.doit().to_sympy()==expected.doit()
    assert actual.tomatrix().to_sympy()==expected.tomatrix()
    assert actual.tolist()==actual.doit().tolist()
    assert tuple(value.to_sympy() for value in actual.shape)==expected.shape
    assert len(actual)==len(expected) and actual.rank()==expected.rank()
    assert actual.to_sympy()==expected
    assert pickle.loads(pickle.dumps(actual))==actual
    assert actual.func(*actual.args)==actual


def test_symbolic_dimensions_substitution_and_bound_index_capture():
    i,j,n,a=sx.symbols('i j n a')
    value=sx.ArrayComprehension(i+a*j,(i,1,n),(j,1,2))
    assert value.doit() is value
    assert value.free_symbols=={a,n}
    assert value.bound_symbols==value.variables==[i,j]
    with pytest.raises(ValueError):len(value)
    with pytest.raises(ValueError):value.tolist()
    with pytest.raises(ValueError):value.tomatrix()
    assert value.subs({n:2,a:3}).tolist()==[[4,7],[5,8]]
    assert value.subs(i,100)==value
    changed=value.subs(a,i)
    assert changed.variables[0]!=i and changed.free_symbols=={i,n}
    assert changed.subs(n,1).tolist()==[[i+1,2*i+1]]
    assert value.xreplace({n:sx.Integer(1),a:sx.Integer(2)}).tolist()==[[3,5]]


@pytest.mark.parametrize('function',[lambda:7,lambda a:a+1,lambda a:a*a])
def test_lambda_maps_use_the_product_of_index_values(function):
    from sympy.tensor.array import ArrayComprehensionMap
    i,j=sx.symbols('i j');I,J=sp.symbols('i j')
    actual=sx.ArrayComprehensionMap(function,(i,1,3),(j,2,4))
    expected=ArrayComprehensionMap(function,(I,1,3),(J,2,4))
    assert actual.doit().to_sympy()==expected.doit()
    assert actual.to_sympy().doit()==expected.doit()


@pytest.mark.parametrize('values,kind',[
    ([1,2],sx.NumberKind),([0,0],sx.NumberKind),([sx.true,sx.false],sx.BooleanKind),
    ([sx.true,1],sx.UndefinedKind),([],sx.UndefinedKind),
])
def test_array_kind_is_interned_and_matches_the_elements(values,kind):
    for constructor in [sx.Array,sx.MutableDenseNDimArray,sx.ImmutableSparseNDimArray]:
        value=constructor(values)
        assert value.kind is sx.ArrayKind(kind)
        assert pickle.loads(pickle.dumps(value.kind)) is value.kind
        # SymPy 1.14 sparse arrays still inherit Basic's UndefinedKind;
        # native sparse arrays expose the same element metadata as dense ones.
        if constructor is not sx.ImmutableSparseNDimArray:
            assert value.kind.to_sympy()==value.to_sympy().kind
    from symox.tensor.array.ndim_array import ArrayKind
    assert ArrayKind is sx.ArrayKind


def test_invalid_limits_and_materialization_bounds():
    i,n=sx.symbols('i n')
    for limits in [(i,0),(i,2,1),(i,0,i),(i+1,0,2),(i,0,sx.oo)]:
        with pytest.raises((TypeError,ValueError)):
            sx.ArrayComprehension(i,limits)
    with pytest.raises(ValueError):sx.ArrayComprehension(i,(i,0,2),(i,0,2))
    with pytest.raises(ValueError):sx.ArrayComprehension(i,(i,0,100000)).doit()
    with pytest.raises(ValueError):sx.ArrayComprehensionMap(lambda a,b:a+b,(i,0,2))
    with pytest.raises(ValueError):sx.ArrayComprehensionMap(sx.sin,(i,0,2))
    with pytest.raises(ValueError):sx.ArrayComprehension(i,(i,0,2)).tomatrix()


def test_comprehensions_do_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,name,*args):
        if name=='sympy' or name.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
i,n=symbols('i n')
a=ArrayComprehension(i*i,(i,1,n)).subs(n,3).doit()
assert a.tolist()==[1,4,9] and a.kind is ArrayKind(NumberKind)
assert ArrayComprehensionMap(lambda a:a+1,(i,1,3)).tolist()==[2,3,4]
'''
    subprocess.run([sys.executable,'-c',source],check=True)
