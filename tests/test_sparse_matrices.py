"""Sparse storage, native sparse algebra and dense exact-solver delegation."""
import pickle
import random
import subprocess
import sys

import numpy as np
import pytest
import sympy as sp
import symox as sx
from symox._sparse import SparseMatrix, MutableSparseMatrix, ImmutableSparseMatrix


def equal(actual, expected):
    ref = actual.to_sympy()
    assert ref.shape == expected.shape
    assert all(sp.simplify(a-b) == 0 for a,b in zip(ref,expected))


@pytest.mark.parametrize('constructor', [MutableSparseMatrix, ImmutableSparseMatrix])
def test_exact_matrix_functions_and_svd_preserve_sparse_storage_and_mutability(constructor):
    a = constructor([[1, 1], [0, 1]])
    x = sx.Symbol('x')
    exponential = a.exp()
    resolvent = a.analytic_func(1/(x+2), x)
    for result in (exponential, resolvent):
        assert isinstance(result, constructor)
    equal(exponential, sp.Matrix([[sp.E, sp.E], [0, sp.E]]))
    equal(resolvent, sp.Matrix([[sp.Rational(1, 3), -sp.Rational(1, 9)], [0, sp.Rational(1, 3)]]))
    rectangular = constructor([[1, sx.I, 2], [sx.I, 1, 0]])
    left, diagonal, right = rectangular.singular_value_decomposition()
    assert all(isinstance(factor, constructor) for factor in (left, diagonal, right))
    equal(left.H*left, sp.eye(2))
    equal(right.H*right, sp.eye(2))
    equal(left*diagonal*right.H, rectangular.to_sympy())


@pytest.mark.parametrize('method', ['exp', 'analytic_func', 'singular_value_decomposition'])
def test_exact_function_limits_precede_sparse_densification(method, monkeypatch):
    a = SparseMatrix(10**8, 10**8, {(0, 0): 1})
    def forbidden(*args, **kwargs):
        raise AssertionError('Oversized sparse input must be rejected before densification')
    monkeypatch.setattr(MutableSparseMatrix, 'as_dense', forbidden)
    with pytest.raises(NotImplementedError, match='before sparse densification'):
        getattr(a, method)()


@pytest.mark.parametrize('args', [
    (), ([[1,0],[0,2]],), ([1,0,2],), ([[1],[0,2]],),
    (2,3,{(0,1):2,(1,2):3}),
    (None,{(3,4):2}),
    (2,2,[1,0,0,2]),
    (2,2,lambda i,j: i+j),
    (0,3,{}), (3,0,{}),
])
def test_sparse_construction(args):
    equal(SparseMatrix(*args),sp.SparseMatrix(*args))


def test_sparse_sympy_numpy_dense_conversion():
    assert sx.SparseMatrix is SparseMatrix
    assert sx.MutableSparseMatrix is MutableSparseMatrix
    assert sx.ImmutableSparseMatrix is ImmutableSparseMatrix
    for source in [sp.SparseMatrix([[1,0],[0,2]]),sx.Matrix([[1,0],[0,2]]),np.array([[1,0],[0,2]])]:
        matrix=SparseMatrix(source)
        equal(matrix,sp.Matrix([[1,0],[0,2]]))
        assert matrix.nnz() == 2
        equal(matrix.as_dense(),sp.Matrix([[1,0],[0,2]]))
        equal(sx.Matrix(matrix),sp.Matrix([[1,0],[0,2]]))
        assert np.array(matrix,dtype=object).shape == (2,2)
    assert SparseMatrix(np.zeros((0,3))).shape == (0,3)
    assert isinstance(sp.sympify(ImmutableSparseMatrix.eye(2)),sp.ImmutableSparseMatrix)
    block=SparseMatrix(4,4,{(1,1):[[1,2],[3,4]]})
    equal(block,sp.Matrix([[0,0,0,0],[0,1,2,0],[0,3,4,0],[0,0,0,0]]))


def test_large_sparse_algebra_never_materializes_dense(monkeypatch):
    size=10**8
    a=SparseMatrix(size,size,{(1,2):3,(7,9):4})
    b=SparseMatrix(size,size,{(2,6):5,(9,0):-2})
    def forbidden(*args,**kwargs): raise AssertionError('Unexpected dense expansion')
    monkeypatch.setattr(MutableSparseMatrix,'as_dense',forbidden)
    result=a*b
    assert result.shape == (size,size)
    assert result.todok() == {(1,6):15,(7,0):-8}
    assert (a+a).todok() == {(1,2):6,(7,9):8}
    assert (a-a).nnz() == 0
    assert a.T.todok() == {(2,1):3,(9,7):4}
    assert (2*a).nnz() == 2
    assert a[1:8,2:10].shape == (7,8)


def test_random_sparse_products_and_cancellations():
    rng=random.Random(813)
    for rows,inner,cols in [(3,4,2),(0,3,2),(3,0,4),(2,3,0),(6,6,6)]:
        data={(i,j):rng.randint(-3,3) for i in range(rows) for j in range(inner) if rng.random()<.3}
        rhs={(i,j):rng.randint(-3,3) for i in range(inner) for j in range(cols) if rng.random()<.3}
        a,b=SparseMatrix(rows,inner,data),SparseMatrix(inner,cols,rhs)
        equal(a*b,a.to_sympy()*b.to_sympy())
        assert (a-a).nnz() == 0
        assert all(v != 0 for v in (a*b).todok().values())


def test_symbolic_scalar_arithmetic_and_substitution():
    x,y=sx.symbols('x y')
    a=SparseMatrix([[x,0],[0,y]])
    equal(a*2,sp.diag(2*sp.Symbol('x'),2*sp.Symbol('y')))
    equal(x*a,a.to_sympy()*sp.Symbol('x'))
    equal(a/2,a.to_sympy()/2)
    assert a.subs({x:0}).nnz() == 1
    assert a.diff(x).todok() == {(0,0):1}
    assert a.xreplace({x:y}).todok() == {(0,0):y,(1,1):y}
    assert a.applyfunc(lambda v:v+1).nnz() == 4
    assert a.free_symbols == {x,y}


@pytest.mark.parametrize('immutable',[False,True])
def test_indexing_slices_and_reshape(immutable):
    cls=ImmutableSparseMatrix if immutable else SparseMatrix
    a=cls([[1,0,2],[0,3,0]])
    assert a[-1] == 0 and a[-1,-2] == 3
    assert a[:] == [1,0,2,0,3,0]
    equal(a[::-1,::-1],a.to_sympy()[::-1,::-1])
    equal(a[:,1],a.to_sympy()[:,1])
    equal(a.reshape(3,2),a.to_sympy().reshape(3,2))
    equal(a.extract([1,0,1],[2,0]),a.to_sympy().extract([1,0,1],[2,0]))
    assert a.row_list() == [(0,0,1),(0,2,2),(1,1,3)]
    assert a.col_list() == [(0,0,1),(1,1,3),(0,2,2)]


def test_mutation_copies_and_sparse_zero_removal():
    a=SparseMatrix.eye(3);b=a.copy()
    a[0,0]=0; assert a.nnz()==2 and b.nnz()==3
    a[0,:]=SparseMatrix([[1,2,0]])
    a[:,2]=SparseMatrix([3,4,5])
    equal(a,sp.Matrix([[1,2,3],[0,1,4],[0,0,5]]))
    a.row_swap(0,2); a.col_swap(0,1); a.row_del(1);a.col_del(2)
    equal(a,sp.Matrix([[0,0],[2,1]]))
    with pytest.raises(TypeError): hash(a)


def test_immutable_storage_and_pickle():
    x=sx.Symbol('x',positive=True)
    a=ImmutableSparseMatrix(3,4,{(1,2):x})
    assert hash(a)==hash(ImmutableSparseMatrix(a))
    for action in [lambda:a.__setitem__((0,0),1),lambda:a.row_swap(0,1),lambda:setattr(a,'_rows',7),lambda:a._data.__setitem__((0,0),sx.Integer(1)),lambda:a.__init__([[2]]),lambda:delattr(a,'_rows')]:
        with pytest.raises((TypeError,AttributeError)): action()
    restored=pickle.loads(pickle.dumps(a))
    assert restored==a and hash(restored)==hash(a)
    assert restored[1,2].is_positive is True
    assert isinstance(a.T,ImmutableSparseMatrix)
    assert isinstance(a+SparseMatrix.zeros(3,4),ImmutableSparseMatrix)
    assert isinstance(a.as_mutable(),MutableSparseMatrix)


def test_sparse_dense_mixing_and_exact_solvers():
    a=SparseMatrix([[2,1],[0,3]])
    dense=sx.eye(2)
    for value in [a+dense,dense+a]:
        assert isinstance(value,sx.Matrix)
        equal(value,a.to_sympy()+sp.eye(2))
    for value in [a*dense,dense*a]:
        assert isinstance(value,sx.Matrix)
        equal(value,a.to_sympy())
    assert isinstance(a.inv(),SparseMatrix)
    equal(a.inv(),a.to_sympy().inv())
    equal(a.LUsolve(SparseMatrix([1,2])),a.to_sympy().inv()*sp.Matrix([1,2]))
    c,f=a.rank_decomposition()
    equal(c*f,a.to_sympy())
    immutable=ImmutableSparseMatrix(a)
    assert isinstance(immutable.inv(),ImmutableSparseMatrix)
    assert isinstance(immutable+dense,sx.ImmutableMatrix)
    assert isinstance(dense.multiply(immutable),sx.ImmutableMatrix)
    assert isinstance(dense.multiply_elementwise(immutable),sx.ImmutableMatrix)


def test_rectangular_sparse_dense_operand_directions():
    sparse=SparseMatrix([[1,0,2],[0,3,0]])
    addend=sx.Matrix([[0,1,0],[2,0,3]])
    for value in [sparse+addend,addend+sparse]:
        equal(value,sparse.to_sympy()+addend.to_sympy())
    equal(addend-sparse,addend.to_sympy()-sparse.to_sympy())
    equal(sparse-addend,sparse.to_sympy()-addend.to_sympy())
    factor=sx.Matrix([[1,0],[0,2],[3,0]])
    equal(sparse*factor,sparse.to_sympy()*factor.to_sympy())
    equal(factor*sparse,factor.to_sympy()*sparse.to_sympy())


def test_sparse_joins_powers_and_numpy():
    a=SparseMatrix([[1,2],[0,1]])
    equal(a**5,a.to_sympy()**5)
    equal(a**-2,a.to_sympy()**-2)
    equal(a**0,sp.eye(2))
    equal(SparseMatrix.hstack(a,a),a.to_sympy().row_join(a.to_sympy()))
    equal(SparseMatrix.vstack(a,a),a.to_sympy().col_join(a.to_sympy()))
    equal(SparseMatrix.diag(a,3),sp.diag(a.to_sympy(),3))
    array=np.array(a,dtype=object)
    assert (np.dot(array,array)==np.array((a*a).tolist(),dtype=object)).all()


def test_invalid_sparse_shapes_indices_and_inputs():
    with pytest.raises(ValueError): SparseMatrix(-1,2,{})
    with pytest.raises(ValueError): SparseMatrix(2,2,{(2,0):1})
    with pytest.raises(ValueError): SparseMatrix(2,2,{(-1,0):1})
    with pytest.raises(TypeError): SparseMatrix(2,2,{0:1})
    with pytest.raises(ValueError): SparseMatrix(2,2,[1])
    with pytest.raises(sx.ShapeError): SparseMatrix.zeros(2)+SparseMatrix.zeros(3)
    with pytest.raises(sx.ShapeError): SparseMatrix.zeros(2,3)*SparseMatrix.zeros(2,3)
    with pytest.raises(IndexError): SparseMatrix.eye(2)[4]
    with pytest.raises(sx.NonSquareMatrixError): SparseMatrix.zeros(2,3)**2


def test_sparse_and_exact_spectral_algorithms_need_no_sympy_runtime():
    source = '''
import sys
import symox as s
a=s.SparseMatrix([[1,2],[0,3]])
assert (a*a).todok() == {(0,0):1,(0,1):8,(1,1):9}
p,j=a.jordan_form()
assert a*p == p*j
assert a.pinv(method='ED') == a.inv()
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable,'-c',source],check=True,capture_output=True,text=True)


@pytest.mark.parametrize('data', [[[1,0],[0,2]], [[0,0],[0,0]], [[sx.Symbol('x'),0]]])
def test_dense_sparse_equality_and_shared_immutable_hash(data):
    sparse=ImmutableSparseMatrix(data)
    dense=sx.ImmutableMatrix(data)
    assert sparse == dense and dense == sparse
    assert not sparse != dense and not dense != sparse
    assert hash(sparse) == hash(dense)
    assert len({sparse,dense}) == 1
    assert {dense:'value'}[sparse] == 'value'
    assert {sparse:'value'}[dense] == 'value'
    assert dense == SparseMatrix(data)
    assert sparse == sx.Matrix(data)


def test_large_immutable_sparse_hash_does_not_materialize_zeros(monkeypatch):
    a=ImmutableSparseMatrix(10**8,10**8,{(1,3):2,(6,8):7})
    def forbidden(*args,**kwargs): raise AssertionError('Unexpected dense expansion')
    monkeypatch.setattr(ImmutableSparseMatrix,'as_dense',forbidden)
    assert hash(a) == hash(ImmutableSparseMatrix(a))
    assert hash(a) != hash(ImmutableSparseMatrix(10**8,10**8,{(1,3):2}))


@pytest.mark.parametrize('shape',[(0,0),(0,3),(3,0)])
def test_empty_dense_sparse_hash(shape):
    sparse=ImmutableSparseMatrix.zeros(*shape)
    dense=sx.ImmutableMatrix.zeros(*shape)
    assert sparse == dense and dense == sparse
    assert hash(sparse) == hash(dense)


@pytest.mark.parametrize('value',[sx.oo,sx.zoo,sx.nan])
def test_sparse_nonfinite_products_include_implicit_zeros(value):
    a=SparseMatrix(2,2,{(0,1):2})
    b=SparseMatrix(2,3,{(0,1):value,(1,2):3})
    result=a*b
    assert result == SparseMatrix(a.to_sympy()*b.to_sympy())
    result=b.T*a.T
    assert result == SparseMatrix(b.T.to_sympy()*a.T.to_sympy())


def test_held_nonfinite_arithmetic_preserves_zero_product_semantics():
    x=sx.Symbol('x')
    for value in [sx.Mul(sx.oo,x,evaluate=False),sx.Add(sx.oo,x,evaluate=False)]:
        special=SparseMatrix([[value]])
        assert (SparseMatrix.zeros(1)*special)[0] == sx.nan
        assert (special*SparseMatrix.zeros(1))[0] == sx.nan
