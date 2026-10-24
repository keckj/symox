"""Native diagonal maps: exact shapes, entries and interoperability.

Contract differences justified by independent shape/entry checks:
* DiagMatrix accepts proven row/column vectors and the empty 0x0 matrix;
  reject nonvectors. SymPy 1.14 accepts these but advertises inconsistent shapes.
* DiagonalOf(A[n,m]) has shape (Min(n,m),1), whereas SymPy currently gives
  (None,1) when dimensions are incomparable. Min(integer,integer) must be integer.
* Preserve declared empty shapes in as_explicit; SymPy sometimes loses 0x1/0xN.
* Symbolic diagonal masking may use native Piecewise(Eq(i,j)) rather than
  introducing an undefined function masquerading as KroneckerDelta.
"""
import pickle

import pytest
import sympy as sp
import symox as sx


CLASSES=('DiagMatrix','DiagonalOf','DiagonalMatrix')


def constructor(module,name):
    return getattr(module,name)


def explicit_same(actual,expected):
    actual=actual.as_explicit().to_sympy()
    assert actual.shape==expected.shape
    assert all(sp.simplify(a-b)==0 for a,b in zip(actual,expected))


@pytest.mark.parametrize('shape',[(3,1),(1,3),(1,1)])
@pytest.mark.parametrize('name',CLASSES)
def test_concrete_interfaces_match_sympy(name,shape):
    A=sx.MatrixSymbol('A',*shape);a=sp.MatrixSymbol('A',*shape)
    value=constructor(sx,name)(A);oracle=constructor(sp,name)(a)
    assert value.shape==oracle.shape
    assert value.args==(A,)
    assert value.func(*value.args)==value
    assert isinstance(value,constructor(sx,name))
    assert value.free_symbols=={A}
    explicit_same(value,oracle.as_explicit())
    assert value.to_sympy()==oracle
    assert sx.MatrixExpr(oracle)==value
    assert sx.sympify(sx.srepr(value))==value


@pytest.mark.parametrize('shape,declared',[
    ((0,1),(0,0)),((1,0),(0,0)),((0,0),(0,0)),
])
def test_empty_vectors_diagonalize_consistently(shape,declared):
    value=sx.DiagMatrix(sx.MatrixSymbol('v',*shape))
    assert value.shape==declared
    assert value.as_explicit().shape==declared
    assert len(value.as_explicit())==0


@pytest.mark.parametrize('shape',[(0,0),(0,1),(1,0),(0,3),(3,0)])
def test_rectangular_diagonal_maps_preserve_empty_shapes(shape):
    A=sx.MatrixSymbol('A',*shape)
    diagonal=sx.DiagonalMatrix(A)
    column=sx.DiagonalOf(A)
    assert diagonal.shape==shape
    assert diagonal.as_explicit().shape==shape
    assert column.shape==(0,1)
    assert column.as_explicit().shape==(0,1)
    assert column.T.as_explicit().shape==(1,0)


@pytest.mark.parametrize('shape',[(2,3),(3,2),(0,3),(3,0)])
def test_diagmatrix_rejects_nonvector_inconsistent_shape_cases(shape):
    with pytest.raises(ValueError):sx.DiagMatrix(sx.MatrixSymbol('M',*shape))


def test_symbolic_dimensions_minimum_is_integer_and_not_an_unknown_shape():
    n,m=sx.symbols('n m',integer=True,nonnegative=True)
    A=sx.MatrixSymbol('A',n,m)
    value=sx.DiagonalOf(A)
    assert value.shape==(sx.Min(n,m),1)
    assert value.rows.is_integer is True
    assert value.rows.is_nonnegative is True
    assert value.free_symbols=={A}
    assert value.xreplace({n:2,m:3}).shape==(2,1)
    assert value.xreplace({n:0,m:3}).as_explicit().shape==(0,1)
    for vector in [sx.MatrixSymbol('v',n,1),sx.MatrixSymbol('w',1,n)]:
        assert sx.DiagMatrix(vector).shape==(n,n)
    with pytest.raises(ValueError):sx.DiagMatrix(A)


@pytest.mark.parametrize('name',CLASSES)
def test_complex_values_transpose_adjoint_and_scalar_derivative(name):
    x=sx.Symbol('x',real=True)
    # A true vector works for all three constructors.
    A=sx.MatrixExpr(sx.ImmutableMatrix([[x+sx.I],[2*x-sx.I],[3]]))
    a=sp.ImmutableMatrix([[sp.Symbol('x',real=True)+sp.I],[2*sp.Symbol('x',real=True)-sp.I],[3]])
    value=constructor(sx,name)(A)
    # SymPy's DiagonalOf/DiagonalMatrix on explicit matrices currently fail
    # when Python int dimensions reach their .as_explicit implementation.
    if name=='DiagMatrix':oracle=sp.diag(*list(a))
    elif name=='DiagonalOf':oracle=sp.ImmutableMatrix(min(a.rows,a.cols),1,lambda i,j:a[i,i])
    else:oracle=sp.ImmutableMatrix(a.rows,a.cols,lambda i,j:a[i,j] if i==j else 0)
    explicit_same(value,oracle)
    explicit_same(value.T,oracle.T)
    explicit_same(value.H,oracle.H)
    explicit_same(value.conjugate(),oracle.conjugate())
    explicit_same(value.diff(x),oracle.diff(sp.Symbol('x',real=True)))
    assert value.diff(sx.Symbol('unrelated'))==sx.ZeroMatrix(*value.shape)


def test_rectangular_diagonal_extraction_and_projection_are_distinct():
    A=sx.MatrixSymbol('A',2,3)
    column=sx.DiagonalOf(A)
    projection=sx.DiagonalMatrix(A)
    assert column.shape==(2,1)
    assert projection.shape==(2,3)
    assert column[1,0]==A[1,1]
    assert projection[1,2]==0
    assert projection[1,1]==A[1,1]
    explicit_same(projection.T,sp.DiagonalMatrix(sp.MatrixSymbol('A',2,3)).as_explicit().T)
    for value in [column,projection]:
        for key in [(value.rows,0),(0,value.cols)]:
            with pytest.raises(IndexError):value[key]


def test_symbolic_entry_masking_matches_dense_values():
    i,j=sx.symbols('i j',integer=True)
    v=sx.MatrixSymbol('v',3,1)
    A=sx.MatrixSymbol('A',3,4)
    diagonal=sx.DiagMatrix(v)[i,j]
    projection=sx.DiagonalMatrix(A)[i,j]
    # Check the predicate mathematically rather than require one printer form.
    for row in range(3):
        for col in range(3):
            assert diagonal.xreplace({i:row,j:col})==(v[row,0] if row==col else 0)
            assert projection.xreplace({i:row,j:col})==(A[row,col] if row==col else 0)
    assert sx.DiagonalOf(A)[i,0]==A[i,i]


def test_one_by_one_vector_orientation_and_negative_indices():
    v=sx.MatrixExpr(sx.ImmutableMatrix([[7]]))
    for name in CLASSES:
        value=constructor(sx,name)(v)
        assert value.shape==(1,1)
        assert value.as_explicit()==sx.ImmutableMatrix([[7]])
        assert value[-1,-1]==7
    row=sx.MatrixExpr(sx.ImmutableMatrix([[2,3,4]]))
    col=sx.MatrixExpr(sx.ImmutableMatrix([[2],[3],[4]]))
    assert sx.DiagMatrix(row).as_explicit()==sx.DiagMatrix(col).as_explicit()
    assert sx.DiagMatrix(row)[-1,-1]==4


@pytest.mark.parametrize('protocol',range(pickle.HIGHEST_PROTOCOL+1))
@pytest.mark.parametrize('name',CLASSES)
def test_pickle_and_dimension_substitution(name,protocol):
    n=sx.Symbol('n',integer=True,positive=True)
    A=sx.MatrixSymbol('A',n,1)
    value=constructor(sx,name)(A)
    restored=pickle.loads(pickle.dumps(value,protocol=protocol))
    assert restored==value and hash(restored)==hash(value)
    assert restored.args==value.args
    assert restored.xreplace({n:3}).as_explicit().shape==constructor(sx,name)(sx.MatrixSymbol('A',3,1)).shape


def test_codegen_numpy_source_without_optional_tensor_packages():
    import numpy as np
    v=sx.MatrixSymbol('v',3,1)
    numeric=np.array([[1.],[2.],[3.]])
    expression=sx.DiagMatrix(v).as_explicit()
    source='numpy.array(['+','.join('['+','.join(sx.numpycode(expression[i,j]) for j in range(3))+']' for i in range(3))+'])'
    np.testing.assert_array_equal(eval(source,{'numpy':np,'v':numeric}),np.diag(numeric[:,0]))


def test_rectangular_projection_normalizes_negative_indices_before_masking():
    A=sx.MatrixSymbol('A',2,3)
    value=sx.DiagonalMatrix(A)
    assert value[-1,-1]==0
    assert value[-1,-2]==A[1,1]
    assert sx.DiagonalOf(A)[-1,0]==A[1,1]


def test_min_max_integer_assumption_is_conservative():
    i,j=sx.symbols('i j',integer=True)
    x=sx.Symbol('x',real=True)
    assert sx.Min(i,j).is_integer is True
    assert sx.Max(i,j).is_integer is True
    assert sx.Min(i,x).is_integer is None
    assert sx.Max(i,x).is_integer is None


def test_diagonal_import_path_and_invalid_vector_replacement():
    from symox.matrices.expressions.diagonal import DiagMatrix,DiagonalOf,DiagonalMatrix
    assert DiagMatrix is sx.DiagMatrix and DiagonalOf is sx.DiagonalOf and DiagonalMatrix is sx.DiagonalMatrix
    vector=sx.MatrixSymbol('v',3,1)
    with pytest.raises(ValueError):sx.DiagMatrix(vector).xreplace({vector:sx.MatrixSymbol('B',2,3)})


def test_diagonal_diff_cache_preserves_alpha_equal_binder_names():
    x,y,z,t=sx.symbols('x y z t')
    f=sx.Function('f')
    # A nonlinear replacement cannot rename the inner derivative variable;
    # this Subs genuinely remains held during differentiation.
    first=sx.Subs(sx.Derivative(f(x),x,evaluate=False),x,t*t+1)
    second=sx.Subs(sx.Derivative(f(y),y,evaluate=False),y,t*t+1)
    assert first==second and first.args!=second.args
    left=sx.DiagMatrix(sx.MatrixExpr(sx.ImmutableMatrix([[first]])))
    right=sx.DiagMatrix(sx.MatrixExpr(sx.ImmutableMatrix([[second]])))
    result=sx.BlockMatrix([[left,right]]).diff(t).as_explicit()
    # SymPy's Subs._eval_derivative form: 2*t*Subs(Derivative(f(x), (x, 2)), x, t**2 + 1).
    held=lambda entry:[a for a in entry.args if a.func==sx.Subs][0]
    assert held(result[0,0]).args[1].args==(x,)
    assert held(result[0,1]).args[1].args==(y,)


def test_diagonal_native_runtime_without_sympy():
    import subprocess
    import sys
    source="""
import sys
sys.modules['sympy']=None
import symox as s
v=s.MatrixSymbol('v',3,1)
assert s.DiagMatrix(v)[1,1]==v[1,0]
assert s.DiagMatrix(v)[1,0]==0
A=s.MatrixSymbol('A',2,3)
assert s.DiagonalOf(A).shape==(2,1)
assert s.DiagonalMatrix(A).as_explicit().shape==(2,3)
"""
    result=subprocess.run([sys.executable,'-c',source],capture_output=True,text=True,timeout=10)
    assert result.returncode==0,result.stderr
