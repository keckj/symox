"""Concrete array protocols, storage ownership and sparse resource bounds."""
import itertools
import pickle
import subprocess
import sys

import numpy as np
import pytest
import sympy as sp
import symox as sx
from symox.tensor.array import (Array,NDimArray,DenseNDimArray,SparseNDimArray,
    ImmutableNDimArray,MutableNDimArray,ImmutableDenseNDimArray,MutableDenseNDimArray,
    ImmutableSparseNDimArray,MutableSparseNDimArray)


KINDS=[ImmutableDenseNDimArray,MutableDenseNDimArray,ImmutableSparseNDimArray,MutableSparseNDimArray]
SHAPES=[(),(1,),(3,),(2,3),(2,2,2),(1,2,1,2)]


def reference(array):
    return getattr(sp,type(array).__name__)([value.to_sympy() for value in array._flat_values()],array.shape)


@pytest.mark.parametrize('kind',KINDS)
@pytest.mark.parametrize('shape',SHAPES)
@pytest.mark.parametrize('source',['flat','nested','numpy','array'])
def test_constructor_shapes_and_ownership(kind,shape,source):
    size=int(np.prod(shape,dtype=np.int64)) if shape else 1
    values=[sx.Integer(i) for i in range(size)]
    if source=='flat':
        array=kind(values,shape)
    elif source=='nested':
        array=kind(Array(values,shape).tolist())
    elif source=='numpy':
        array=kind(np.array(list(range(size))).reshape(shape))
    else:
        array=kind(MutableDenseNDimArray(values,shape))
    assert array.shape==shape
    assert array.rank()==len(shape) and len(array)==size
    expected=reference(array)
    assert list(array._flat_values())==values
    assert array.to_sympy()==expected
    assert array.func(*array.args)==array
    assert pickle.loads(pickle.dumps(array))==array
    if shape:
        assert len(list(array))==shape[0]
    assert array.copy()==array
    assert array.as_immutable()==array and array.as_mutable()==array
    assert isinstance(array.as_immutable(),ImmutableNDimArray)
    assert isinstance(array.as_mutable(),MutableNDimArray)
    if not array._mutable:
        assert hash(array)==hash(array.as_immutable())
        with pytest.raises(TypeError):
            array[(0,)*array.rank()]=sx.Integer(9)
    else:
        with pytest.raises(TypeError):
            hash(array)
        copy=array.copy()
        copy[(0,)*array.rank()]=sx.Integer(9)
        assert array[(0,)*array.rank()]==0


@pytest.mark.parametrize('kind',KINDS)
@pytest.mark.parametrize('key',[
    (),0,-1,(1,), (slice(None),1),(slice(None,None,-1),slice(0,3,2)),
    (1,slice(None,None,-1)),(slice(0,0),slice(None)),
])
def test_indices_and_slices_match_reference(kind,key):
    array=kind(list(range(6)),(2,3))
    expected=reference(array)[key]
    if key==(slice(None,None,-1),slice(0,3,2)):
        # SymPy 1.14 truncates a positive stride using floor division. Python
        # slicing includes both column 0 and column 2 in this example.
        expected=sp.Array([[3,5],[0,2]])
    actual=array[key]
    assert actual.to_sympy()==expected
    if isinstance(actual,NDimArray):
        assert actual.shape==tuple(expected.shape)
        assert isinstance(actual,kind)


@pytest.mark.parametrize('kind',KINDS)
def test_nested_tuple_and_mapping_inputs(kind):
    x=sx.Symbol('x')
    value=kind([sx.Tuple(1,x),sx.Tuple(3,4)])
    assert value.shape==(2,2) and value[0,1]==x
    mapping={(0,0):x,3:2}
    array=kind(mapping,(2,2))
    assert mapping=={(0,0):x,3:2}
    assert array.tolist()==[[x,0],[0,2]]
    matrix=sx.Matrix([[x,0],[0,2]])
    assert kind(matrix)==array
    assert kind(sx.SparseMatrix(matrix))==array
    assert kind(sx.MatrixSymbol('A',2,2)).shape==(2,2)


@pytest.mark.parametrize('kind',KINDS)
def test_algebra_and_native_scalar_coefficients(kind):
    x,y=sx.symbols('x y')
    array=kind([[x,2],[0,y]])
    other=kind([[y,1],[3,x]])
    expected,ref_other=reference(array),reference(other)
    for actual,ref in [(array+other,expected+ref_other),(array-other,expected-ref_other),
        (-array,-expected),(2*array,2*expected),(array/2,expected/2)]:
        assert actual.to_sympy()==ref
    assert sx.QQ(2)*array==array*2
    assert sx.QQ_I.convert(sx.I)*array==array*sx.I
    assert array.free_symbols=={x,y}
    assert array.subs(x,y).tolist()==[[y,2],[0,y]]
    assert array.xreplace({x:y}).tolist()==[[y,2],[0,y]]
    assert array.diff(x).tolist()==[[1,0],[0,0]]
    assert not array.diff(x)._mutable
    assert array.applyfunc(lambda v:v+1).tolist()==[[x+1,3],[1,y+1]]
    assert array.T.tolist()==[[x,0],[2,y]]
    assert array.tomatrix()==sx.Matrix([[x,2],[0,y]])
    assert array.reshape(4).tolist()==[x,2,0,y]
    with pytest.raises(ValueError):
        array*other
    with pytest.raises(ValueError):
        array+kind([1,2])


@pytest.mark.parametrize('kind',[MutableDenseNDimArray,MutableSparseNDimArray])
def test_mutation_is_owned_and_failed_slice_assignment_is_atomic(kind):
    source=[[1,2],[3,4]]
    array=kind(source)
    source[0][0]=100
    assert array[0,0]==1
    row=array[0]
    row[0]=200
    assert array[0,0]==1
    array[:,1]=[7,8]
    assert array.tolist()==[[1,7],[3,8]]
    array[0]=[5,6]
    assert array.tolist()==[[5,6],[3,8]]
    previous=array.copy()
    for replacement in [3,[1,2,3],[[1,2]]]:
        with pytest.raises(ValueError):
            array[:,0]=replacement
        assert array==previous
    vector=kind([1,2,3])
    vector[::-1]=vector
    assert vector.tolist()==[3,2,1]
    assert vector.as_immutable()==ImmutableDenseNDimArray([3,2,1])


@pytest.mark.parametrize('kind',KINDS)
def test_empty_axes_and_scalar_tensors_have_consistent_shapes(kind):
    for shape in [(0,),(0,2),(2,0),(2,0,3)]:
        array=kind.zeros(*shape)
        assert array.shape==shape and len(array)==0
        assert array.copy().shape==shape
        assert array.reshape(0).shape==(0,)
        assert pickle.loads(pickle.dumps(array))==array
    assert kind.zeros(2,0).tolist()==[[],[]]
    assert kind.zeros(0,2).tolist()==[]
    scalar=kind([3],())
    assert scalar[()]==3 and scalar.tolist()==3 and list(scalar)==[sx.Integer(3)]
    assert kind.zeros()[()]==0 and kind.ones()[()]==1


def test_sparse_large_shapes_never_densify_in_structural_operations():
    x=sx.Symbol('x')
    array=ImmutableSparseNDimArray({(1,3):x},(2,10**10))
    assert array[1,3]==x and array[0,3]==0
    assert array[1].shape==(10**10,)
    assert array[1][3]==x
    assert array[:,3].tolist()==[0,x]
    assert array[:,::-1][1,10**10-4]==x
    assert array.reshape(10**10,2).shape==(10**10,2)
    assert array.copy()==array
    assert array.free_symbols=={x}
    assert pickle.loads(pickle.dumps(array))==array
    assert len(array._data)==1
    assert (array+array)[1,3]==2*x
    with pytest.raises(ValueError):
        array.tolist()
    with pytest.raises(ValueError):
        array.as_mutable().applyfunc(lambda value:value+1)
    empty=ImmutableDenseNDimArray([], (10**10,0))
    assert empty[0].shape==(0,)
    assert empty[:].shape==empty.shape
    with pytest.raises(ValueError):
        empty.tolist()
    matrix=sx.SparseMatrix(10**10,2,{(1,1):x})
    converted=ImmutableSparseNDimArray(matrix, (2,10**10))
    assert converted[0,3]==x
    with pytest.raises(ValueError):
        Array(matrix)


def test_equal_immutable_storage_kinds_have_equal_hashes():
    dense=ImmutableDenseNDimArray([[1,0],[0,2]])
    sparse=ImmutableSparseNDimArray({0:1,3:2},(2,2))
    assert dense==sparse and hash(dense)==hash(sparse)
    assert len({dense,sparse})==1
    assert dense!=[[1,0],[0,2]]


@pytest.mark.parametrize('kind',KINDS)
def test_symbolic_indices_substitute_valid_negative_and_positive_positions(kind):
    i,j=sx.symbols('i j',integer=True)
    array=kind([[1,0,3],[4,5,6]])
    entry=array[i,j]
    row=array[i,:]
    for r,c in itertools.product(range(-2,2),range(-3,3)):
        assert entry.subs({i:r,j:c})==array[r,c]
        assert row.subs(i,r)[c]==array[r,c]
    generic=sx.Symbol('k')
    assert array[generic,0].subs(generic,-1)==4
    assert array[generic,0].subs(generic,sx.Rational(1,2))==sx.nan


def test_numpy_conversion_and_explicit_copy_semantics():
    from symox.tensor.array import Array
    x=sx.Symbol('x')
    array=Array([[1,x],[3,4]])
    objects=np.asarray(array)
    assert objects.shape==(2,2) and objects.dtype==object
    assert objects[0,1]==x
    objects[0,0]=10
    assert array[0,0]==1
    numeric=np.asarray(Array([[1,2],[3,4]]),dtype=float)
    np.testing.assert_array_equal(numeric,[[1,2],[3,4]])
    assert np.asarray(Array([], (0,3))).shape==(0,3)
    with pytest.raises(ValueError):
        array.__array__(copy=False)


def test_invalid_shapes_indices_and_rectangular_input():
    for shape in [(-1,), (1.5,), (2,2)]:
        with pytest.raises((TypeError,ValueError)):
            Array([1,2],shape)
    for source in [[[1,2],[3]],[1,[2,3]]]:
        with pytest.raises(ValueError):
            Array(source)
    cyclic=[]
    cyclic.append(cyclic)
    with pytest.raises(ValueError):
        Array(cyclic)
    array=Array([[1,2],[3,4]])
    for index in [2,-3,(0,0,0),(0,1.5),(slice(None,None,0),0)]:
        with pytest.raises((ValueError,TypeError)):
            array[index]
    with pytest.raises(ValueError):
        array.reshape(3)
    with pytest.raises(ValueError):
        Array.zeros(100001)
    with pytest.raises(ValueError):
        ImmutableSparseNDimArray({(2,0):1},(2,2))
    with pytest.raises(ValueError):
        Array([1,2]).tomatrix()
    with pytest.raises(ValueError):
        Array([1,2]).transpose()


def test_no_sympy_import():
    source='''
import sys,pickle
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import symbols,QQ
from symox.tensor.array import Array,MutableSparseNDimArray,tensorproduct,tensorcontraction
x,y=symbols('x y')
A=Array([[x,0],[0,y]])
assert tensorcontraction(A,(0,1))==x+y
assert tensorproduct(A,A).shape==(2,2,2,2)
assert (QQ(2)*A)[0,0]==2*x
B=MutableSparseNDimArray({},(2,2)); B[1,0]=x
assert pickle.loads(pickle.dumps(B))==B
'''
    subprocess.run([sys.executable,'-c',source],check=True)
