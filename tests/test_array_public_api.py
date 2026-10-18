"""Public dispatch, finite matrix calculus and shape registration."""
import itertools
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('name', ['Array','NDimArray','DenseNDimArray','SparseNDimArray',
    'MutableNDimArray','ImmutableNDimArray','MutableDenseNDimArray',
    'ImmutableDenseNDimArray','MutableSparseNDimArray','ImmutableSparseNDimArray',
    'tensorproduct','tensorcontraction','tensordiagonal','permutedims','derive_by_array'])
def test_root_exports(name):
    from symox.tensor import array
    assert getattr(sx,name) is getattr(array,name)
    assert name in sx.__all__


@pytest.mark.parametrize('name', ['simplify','expand','factor','trigsimp','powsimp',
    'ratsimp','radsimp','cancel','together','conjugate','N'])
@pytest.mark.parametrize('kind', [sx.Array,sx.MutableDenseNDimArray,sx.ImmutableSparseNDimArray,
    sx.Matrix,sx.ImmutableMatrix,sx.SparseMatrix,sx.ImmutableSparseMatrix])
def test_scalar_operations_dispatch_entrywise(name,kind):
    x=sx.Symbol('x')
    array=kind([sx.sin(x)**2+sx.cos(x)**2,(x*x-1)/(x-1),sx.sqrt(2)])
    operation=getattr(sx,name)
    assert operation(array)==array.applyfunc(operation)
    assert isinstance(operation(array),kind)
    assert operation(x).to_sympy()==getattr(sp,name)(sp.Symbol('x'))


def test_public_gradients_hessians_and_explicit_matrix_variables():
    x,y=sx.symbols('x y');variables=sx.Array([x,y])
    value=x*x*y+sx.sin(y)
    assert sx.diff(value,variables)==value.diff(variables)==sx.Array([2*x*y,x*x+sx.cos(y)])
    hessian=sx.diff(sx.diff(value,variables),variables)
    assert hessian==sx.Array([[2*y,2*x],[2*x,-sx.sin(y)]])
    assert sx.diff(sx.Array([x*x,y*y]),x,2)==sx.Array([2,0])
    matrix=sx.Matrix([[x,y]])
    derivative=matrix.diff(sx.Matrix([x,y]))
    assert derivative.shape==(2,1,1,2)
    assert derivative[0,0,0,0]==1 and derivative[1,0,0,1]==1
    assert sx.Array([x,y]).diff(sx.Matrix([x,y])).shape==(2,1,2)


@pytest.mark.parametrize('operation', ['identity','transpose','square','inverse','trace','det'])
def test_finite_matrix_jacobians_match_explicit_scalar_derivatives(operation):
    A=sx.MatrixSymbol('A',2,2);B=sp.MatrixSymbol('A',2,2)
    native={'identity':lambda A:A,'transpose':lambda A:A.T,'square':lambda A:A*A,
        'inverse':lambda A:A.inv(),'trace':sx.Trace,'det':sx.Determinant}[operation](A)
    reference={'identity':lambda A:A,'transpose':lambda A:A.T,'square':lambda A:A*A,
        'inverse':lambda A:A.inv(),'trace':sp.Trace,'det':sp.Determinant}[operation](B)
    reference=reference.as_explicit() if hasattr(reference,'as_explicit') else reference.doit()
    if operation=='det':
        reference=B.as_explicit().det()
    derivative=native.diff(A)
    assert derivative==sx.diff(native,A)
    if hasattr(reference,'shape'):
        for i,j,k,l in itertools.product(range(2),repeat=4):
            assert sp.cancel(derivative[i,j,k,l].to_sympy()-reference[k,l].diff(B[i,j]))==0
    else:
        for i,j in itertools.product(range(2),repeat=2):
            assert sp.cancel(derivative[i,j].to_sympy()-reference.diff(B[i,j]))==0


def test_derivative_limits_and_options():
    x=sx.Symbol('x');variables=sx.Array([x])
    assert sx.diff(x,variables,evaluate=True)==sx.Array([1])
    for value in [x,sx.Array([x]),sx.Matrix([[x]]),sx.MatrixSymbol('A',1,1)]:
        for options in [dict(evaluate=False),dict(unknown=True)]:
            with pytest.raises(NotImplementedError):
                sx.diff(value,variables,**options)
    n=sx.Symbol('n',integer=True,positive=True)
    A=sx.MatrixSymbol('A',n,n)
    with pytest.raises((ValueError,NotImplementedError,TypeError)):
        A.diff(A)


def test_shape_queries_and_registration():
    from symox.tensor.functions import NoShapeError
    for value in [sx.Array([1,2]),sx.Matrix([[1,2]]),sx.MatrixSymbol('A',2,3)]:
        assert sx.shape(value)==value.shape
    with pytest.raises(NoShapeError):
        sx.shape(sx.Symbol('x'))
    class Wrapper:
        def __init__(self,value):self.value=value
    @sx.shape.register(Wrapper)
    def wrapper_shape(value):return sx.shape(value.value)
    assert sx.shape(Wrapper(sx.Array([1,2])))==(2,)


def test_no_sympy_import_for_root_array_api():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x,y=symbols('x y')
assert diff(x*x*y,Array([x,y])).tolist()==[2*x*y,x*x]
assert simplify(Array([sin(x)**2+cos(x)**2])).tolist()==[1]
A=MatrixSymbol('A',2,2)
assert Trace(A).diff(A).tolist()==[[1,0],[0,1]]
assert A.diff(A)[0,1,0,1]==1
assert shape(tensorproduct(Array([1,2]),Array([3,4])))==(2,2)
'''
    subprocess.run([sys.executable,'-c',source],check=True)
