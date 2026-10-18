"""Tensor identities and differential checks across concrete storage classes."""
import itertools
import math
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from symox.tensor.array import (Array, NDimArray, ImmutableDenseNDimArray,
    MutableDenseNDimArray, ImmutableSparseNDimArray, MutableSparseNDimArray,
    tensorproduct,tensorcontraction,tensordiagonal,permutedims,derive_by_array)


KINDS=[ImmutableDenseNDimArray,MutableDenseNDimArray,ImmutableSparseNDimArray,MutableSparseNDimArray]


def same(actual,expected):
    if isinstance(actual,NDimArray):
        assert actual.shape==tuple(expected.shape)
        assert list(actual.to_sympy().reshape(len(actual)))==list(expected.reshape(len(actual)))
    else:
        assert sp.simplify(actual.to_sympy()-expected)==0


@pytest.mark.parametrize('first_kind',KINDS)
@pytest.mark.parametrize('second_kind',KINDS)
@pytest.mark.parametrize('shapes',[((2,),(3,)),((2,2),(2,)),((2,1),(1,3))])
def test_tensorproducts_match_reference(first_kind,second_kind,shapes):
    x=sx.Symbol('x')
    first=first_kind([x if i==1 else i for i in range(math.prod(shapes[0]))],shapes[0])
    second=second_kind([i+1 for i in range(math.prod(shapes[1]))],shapes[1])
    actual=tensorproduct(first,second)
    expected=sp.tensorproduct(first.to_sympy(),second.to_sympy())
    same(actual,expected)
    assert not actual._mutable
    for index in itertools.product(*(range(d) for d in actual.shape)):
        assert actual[index]==first[index[:first.rank()]]*second[index[first.rank():]]
    same(tensorproduct(2,first,3,second),sp.tensorproduct(2,first.to_sympy(),3,second.to_sympy()))


@pytest.mark.parametrize('kind',KINDS)
@pytest.mark.parametrize('shape,axes',[
    ((2,2),((0,1),)),((2,3,2),((0,2),)),((2,2,2),((0,1,2),)),
    ((2,3,2,3),((0,2),(1,3))),((2,3),((0,),)),
])
def test_contractions_and_diagonals_match_reference(kind,shape,axes):
    import math
    values=list(range(math.prod(shape)))
    array=kind(values,shape)
    same(tensorcontraction(array,*axes),sp.tensorcontraction(array.to_sympy(),*axes))
    if any(len(group)<2 for group in axes):
        with pytest.raises(ValueError):
            tensordiagonal(array,*axes)
    else:
        same(tensordiagonal(array,*axes),sp.tensordiagonal(array.to_sympy(),*axes))
    assert tensorcontraction(array) is array
    assert tensordiagonal(array) is array


@pytest.mark.parametrize('kind',KINDS)
@pytest.mark.parametrize('permutation',list(itertools.permutations(range(3))))
def test_axis_permutations_and_inverse(kind,permutation):
    array=kind(list(range(12)),(2,3,2))
    actual=permutedims(array,permutation)
    same(actual,sp.permutedims(array.to_sympy(),permutation))
    inverse=tuple(permutation.index(i) for i in range(3))
    assert permutedims(actual,inverse)==array
    old='ijk';new=''.join(old[i] for i in permutation)
    assert permutedims(array,index_order_old=old,index_order_new=new)==actual
    from symox.combinatorics import Permutation
    assert permutedims(array,Permutation(permutation))==actual
    assert isinstance(actual,kind)


@pytest.mark.parametrize('kind',KINDS)
def test_gradients_jacobians_and_higher_derivatives(kind):
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    variables=Array([x,y])
    array=kind([x*x+y,sx.sin(x*y)])
    expected=sp.Array([X*X+Y,sp.sin(X*Y)])
    same(derive_by_array(array,variables),sp.derive_by_array(expected,sp.Array([X,Y])))
    same(array.diff(variables),sp.derive_by_array(expected,sp.Array([X,Y])))
    same(derive_by_array(x*x*y,[x,y]),sp.derive_by_array(X*X*Y,[X,Y]))
    same(derive_by_array(array,x),sp.derive_by_array(expected,X))
    same(array.diff(x,2),expected.diff(X,2))
    same(array.diff((x,2)),expected.diff((X,2)))
    hessian=derive_by_array(derive_by_array(x*x*y,[x,y]),[x,y])
    assert hessian.tolist()==[[2*y,2*x],[2*x,0]]
    assert array.diff(x).free_symbols=={x,y}


def test_default_differentiation_and_invalid_variables():
    x,y=sx.symbols('x y')
    assert Array([x,x*x]).diff().tolist()==[1,2*x]
    for array in [Array([x,y]),Array([1,2]),Array([])]:
        with pytest.raises(ValueError):
            array.diff()
    for variable in [0,x*x,sx.sin(x)]:
        with pytest.raises(ValueError):
            derive_by_array(Array([x]),[variable])
        with pytest.raises(ValueError):
            Array([x]).diff(variable)
    f=sx.Function('f')
    assert derive_by_array(f(x)**2,[f(x)])[0]==2*f(x)


def test_finite_matrix_symbols_have_explicit_array_jacobians():
    matrix=sx.MatrixSymbol('A',2,2)
    array=Array(matrix)
    derivative=array.diff(matrix)
    assert derivative.shape==(2,2,2,2)
    for i,j,k,l in itertools.product(range(2),repeat=4):
        assert derivative[i,j,k,l]==int(i==k and j==l)
    x=sx.Symbol('x')
    scalar=sx.trace(matrix)
    gradient=derive_by_array(scalar,matrix)
    assert gradient.tolist()==[[1,0],[0,1]]
    assert Array([scalar]).diff(matrix[0,0]).tolist()==[1]
    determinant=derive_by_array(sx.det(matrix),matrix)
    assert determinant.tolist()==[[matrix[1,1],-matrix[1,0]],[-matrix[0,1],matrix[0,0]]]
    with pytest.raises((TypeError,NotImplementedError)):
        Array(sx.MatrixSymbol('B',x,x))


def test_tensor_contraction_reconstructs_ordered_matrix_multiplication():
    A=sx.MatrixSymbol('A',2,3)
    B=sx.MatrixSymbol('B',3,2)
    product=tensorcontraction(tensorproduct(A,B),(1,2))
    assert product.tomatrix()==(A*B).as_explicit()


def test_sparse_tensor_operations_keep_large_implicit_zero_spaces():
    x,y=sx.symbols('x y')
    a=ImmutableSparseNDimArray({(1,2):x},(10**6,10**6))
    b=ImmutableSparseNDimArray({3:y},(4,))
    product=tensorproduct(a,b)
    assert product.shape==(10**6,10**6,4)
    assert product[1,2,3]==x*y and len(product._data)==1
    assert tensorcontraction(a,(0,1))==0
    diagonal=tensordiagonal(a,(0,1))
    assert diagonal.shape==(10**6,) and len(diagonal._data)==0
    permuted=permutedims(a,(1,0))
    assert permuted[2,1]==x and len(permuted._data)==1
    gradient=derive_by_array(a,[x,y])
    assert gradient.shape==(2,10**6,10**6) and gradient[0,1,2]==1
    assert len(gradient._data)==1
    assert pickle.loads(pickle.dumps(gradient))==gradient


@pytest.mark.parametrize('kind',[ImmutableSparseNDimArray,MutableSparseNDimArray])
def test_nonfinite_tensorproducts_preserve_zero_times_infinity(kind):
    first=kind([0,1]);second=kind([sx.oo,2])
    actual=tensorproduct(first,second)
    assert actual[0,0]==sx.nan and actual[0,1]==0
    assert actual[1,0]==sx.oo and actual[1,1]==2
    with pytest.raises(ValueError):
        tensorproduct(kind({},(10**6,)),sx.oo)


def test_scalar_arguments_and_empty_axes():
    assert tensorproduct()==1 and tensorproduct(2,3)==6
    for kind in KINDS:
        empty=kind.zeros(0,2)
        assert tensorproduct(empty,Array([1,2])).shape==(0,2,2)
        assert tensorcontraction(kind.zeros(0,0),(0,1))==0
        assert tensordiagonal(kind.zeros(0,0),(0,1)).shape==(0,)
        assert permutedims(empty,(1,0)).shape==(2,0)
        assert empty.applyfunc(lambda value:1/int(value)).shape==(0,2)
    scalar=Array([3],())
    assert tensorcontraction(scalar) is scalar


@pytest.mark.parametrize('operation',[tensorcontraction,tensordiagonal])
def test_invalid_axis_groups_are_rejected(operation):
    array=Array(list(range(12)),(2,3,2))
    for groups in [((0,1),),((0,0),),((0,2),(2,)),((3,),),((),),((-1,),)]:
        with pytest.raises(ValueError):
            operation(array,*groups)


def test_invalid_permutations_and_materialization_limits():
    array=Array([[1,2],[3,4]])
    for permutation in [None,[0],[0,0],[0,2],[0,1,2]]:
        with pytest.raises(ValueError):
            permutedims(array,permutation)
    with pytest.raises(ValueError):
        permutedims(array,index_order_old='ii',index_order_new='ii')
    with pytest.raises(ValueError):
        tensorproduct(Array.ones(400),Array.ones(400))


def test_no_sympy_import():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import symbols,MatrixSymbol
from symox.tensor.array import *
x,y=symbols('x y')
A=Array([[x,y],[1,2]])
assert tensorcontraction(A,(0,1))==x+2
assert derive_by_array(x*x*y,[x,y]).tolist()==[2*x*y,x*x]
assert permutedims(A,[1,0]).shape==(2,2)
B=MatrixSymbol('B',2,2)
assert Array(B).diff(B).shape==(2,2,2,2)
'''
    subprocess.run([sys.executable,'-c',source],check=True)
