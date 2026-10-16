"""Python-compatible concrete MatrixSlice semantics.

The normative oracle is list slicing, not SymPy's
MatrixExpr slice implementation (which floors non-divisible lengths and rejects
many valid empty/reverse ranges). Empty results remain 2-D matrices.

The core needs a typed normalized SliceSpec. In particular an omitted reverse
stop is distinct from explicit -1. Reconstruction, composition, codec and srepr
must never normalize an already-normalized -1 sentinel a second time.
"""
import itertools
import operator
import pickle

import pytest
import symox as sx


BOUNDS=(None,-10,-2,-1,0,1,2,10)
STEPS=(-3,-2,-1,1,2,3)
SHAPES=((0,0),(0,4),(4,0),(1,1),(2,3),(3,2),(5,6))


def indices(size,spec):
    return list(range(size))[spec]


def assert_numeric_view(view,parent_shape,rows,cols):
    row_ids=indices(parent_shape[0],rows)
    col_ids=indices(parent_shape[1],cols)
    expected_shape=(len(row_ids),len(col_ids))
    assert view.shape==expected_shape
    result=view.as_explicit()
    assert result.shape==expected_shape
    assert list(result)==[100*i+j for i in row_ids for j in col_ids]


def numeric_parent(rows,cols):
    # Use the MatrixExpr wrapper to exercise the abstract slice DAG rather than
    # the existing concrete DenseMatrix indexing implementation.
    return sx.MatrixExpr(sx.ImmutableMatrix(rows,cols,[100*i+j for i in range(rows) for j in range(cols)]))


@pytest.mark.parametrize('shape',SHAPES)
@pytest.mark.parametrize('axis',[0,1])
def test_cartesian_concrete_bounds_steps_and_empty_shapes(shape,axis):
    parent=numeric_parent(*shape)
    # 5,376 independent combinations across these parameter groups. Avoid
    # multiplying every row case by every column case; two-axis compositions
    # and interleaving are certified separately below.
    for start,stop,step in itertools.product(BOUNDS,BOUNDS,STEPS):
        selected=slice(start,stop,step)
        row,col=(selected,slice(None)) if axis==0 else (slice(None),selected)
        view=parent[row,col]
        assert_numeric_view(view,shape,row,col)


@pytest.mark.parametrize('shape',[(0,4),(4,0),(1,1),(5,6)])
@pytest.mark.parametrize('rows,cols',[
    (slice(None,None,-1),slice(None,None,-1)),
    (slice(None,None,2),slice(None,None,-2)),
    (slice(-100,100,3),slice(100,-100,-2)),
    (slice(2,2),slice(None)),
    (slice(None),slice(3,1)),
    (slice(-1,None,-1),slice(-1,-7,-2)),
])
def test_two_axis_ranges_and_omitted_reverse_stop(shape,rows,cols):
    parent=numeric_parent(*shape)
    view=parent[rows,cols]
    assert_numeric_view(view,shape,rows,cols)
    assert view.func(*view.args)==view
    assert sx.sympify(sx.srepr(view))==view


def test_explicit_negative_one_and_omitted_reverse_stop_are_distinct():
    parent=numeric_parent(5,6)
    reverse=parent[4::-1,:]
    empty=parent[4:-1:-1,:]
    assert reverse.shape==(5,6)
    assert empty.shape==(0,6)
    assert reverse!=empty
    assert_numeric_view(reverse,(5,6),slice(4,None,-1),slice(None))
    for value in [reverse,empty]:
        assert value.func(*value.args)==value
        assert sx.sympify(sx.srepr(value))==value
        assert pickle.loads(pickle.dumps(value))==value


@pytest.mark.parametrize('outer_rows,outer_cols,inner_rows,inner_cols',[
    (slice(None,None,2),slice(None,None,-1),slice(None,None,-1),slice(1,None,2)),
    (slice(None,None,-1),slice(None,None,-2),slice(1,None,2),slice(None,None,-1)),
    (slice(1,5,2),slice(5,0,-2),slice(0,10),slice(-100,100,2)),
    (slice(1,1),slice(None),slice(None,None,-1),slice(None,None,-1)),
    (slice(None),slice(0,0),slice(None,None,-2),slice(None)),
])
def test_nested_slice_composition_preserves_stride_and_empty_shapes(outer_rows,outer_cols,inner_rows,inner_cols):
    parent=numeric_parent(5,6)
    view=parent[outer_rows,outer_cols][inner_rows,inner_cols]
    rows=indices(5,outer_rows)[inner_rows]
    cols=indices(6,outer_cols)[inner_cols]
    assert view.shape==(len(rows),len(cols))
    assert list(view.as_explicit())==[100*i+j for i in rows for j in cols]
    # A composed view should remain attached to the original parent when its
    # normalized arithmetic progression can be represented directly.
    assert view.parent==parent
    assert sx.sympify(sx.srepr(view))==view


@pytest.mark.parametrize('protocol',range(pickle.HIGHEST_PROTOCOL+1))
@pytest.mark.parametrize('row,col',[
    (slice(None,None,-1),slice(None,None,2)),
    (slice(4,-1,-1),slice(None)),
    (slice(None),slice(2,2,-1)),
    (slice(-100,100,3),slice(100,-100,-2)),
])
def test_flat_codec_preserves_normalized_slice_spec_and_omissions(protocol,row,col):
    A=sx.MatrixSymbol('A',5,6)
    value=A[row,col]
    restored=pickle.loads(pickle.dumps(value,protocol=protocol))
    assert restored==value and restored.args==value.args
    assert restored.shape==value.shape
    assert restored.func(*restored.args)==restored
    assert sx.sympify(sx.srepr(restored))==restored


def test_scalar_axis_plus_slice_keeps_a_two_dimensional_view():
    parent=numeric_parent(5,6)
    assert_numeric_view(parent[-1,::-2],(5,6),slice(4,5),slice(None,None,-2))
    assert_numeric_view(parent[::-2,-1],(5,6),slice(None,None,-2),slice(5,6))
    for key in [(5,slice(None)),(-6,slice(None)),(slice(None),6)]:
        with pytest.raises(IndexError):parent[key]


def test_slice_length_uses_ceiling_instead_of_truncation():
    parent=numeric_parent(5,6)
    assert parent[::2,:].shape==(3,6)
    assert parent[-2:-5:-2,:].shape==(2,6)
    assert parent[::-2,:].shape==(3,6)


@pytest.mark.parametrize('row,col',[
    (slice(None,None,-1),slice(None,None,2)),
    (slice(1,1),slice(None)),
])
def test_diagonal_scalar_derivative_and_adjoint_keep_slice_ranges(row,col):
    x=sx.Symbol('x',real=True)
    A=sx.MatrixSymbol('A',5,6)
    view=(x*A)[row,col]
    assert view.diff(x)==A[row,col]
    expected=A[row,col].as_explicit()
    assert view.diff(x).as_explicit()==expected
    assert view.T.as_explicit()==view.as_explicit().T
    assert view.H.as_explicit()==view.as_explicit().H
    assert view.diff(sx.Symbol('y'))==sx.ZeroMatrix(*view.shape)


@pytest.mark.parametrize('spec',[slice(None,None,0),slice(None,None,False)])
def test_zero_step_is_rejected_like_python(spec):
    with pytest.raises(ValueError):numeric_parent(5,6)[spec,:]


@pytest.mark.parametrize('spec',[slice(1.25,None),slice(None,2.5),slice(None,None,1.5)])
def test_known_noninteger_bounds_or_steps_are_rejected(spec):
    with pytest.raises(TypeError):numeric_parent(5,6)[spec,:]


def test_concrete_integer_index_protocol_and_boolean_bounds():
    class Index:
        def __init__(self,value):self.value=value
        def __index__(self):return self.value
    parent=numeric_parent(5,6)
    cases=[slice(True,None,True),slice(True,False,-1),slice(Index(-1),None,Index(-2))]
    for raw in cases:
        spec=slice(*(None if x is None else operator.index(x) for x in (raw.start,raw.stop,raw.step)))
        assert_numeric_view(parent[raw,:],(5,6),spec,slice(None))


def test_forward_sympy_conversion_represents_exact_selected_values():
    # Native -> SymPy needs an explicit policy because SymPy's normal
    # MatrixSlice constructor rejects reverse/empty and computes floor lengths.
    # Require mathematical conversion, without requiring a broken class shape.
    parent=numeric_parent(5,6)
    for rows,cols in [(slice(None,None,-1),slice(None)),(slice(None,None,2),slice(None,None,2)),(slice(1,1),slice(None)),(slice(None),slice(2,2)),(slice(0,0),slice(0,0))]:
        view=parent[rows,cols]
        oracle=view.to_sympy()
        explicit=oracle if hasattr(oracle,'tolist') else oracle.as_explicit()
        assert explicit.shape==view.shape
        assert list(explicit)==[int(v) for v in view.as_explicit()]


def test_symbolic_bounds_steps_and_substitution_match_list_slices():
    n = sx.Symbol('n', integer=True, nonnegative=True)
    k = sx.Symbol('k', integer=True)
    h = sx.Symbol('h', integer=True, positive=True)
    A = sx.MatrixSymbol('A', n, 4)
    for raw in [slice(None, None, h), slice(None, None, -h), slice(k, None, h), slice(None, k, -h)]:
        view = A[raw, :]
        assert view.rows.is_integer is True
        assert view.rows.is_nonnegative is True
        for size, bound, stride in [(0, -9, 2), (1, 8, 3), (5, -2, 2), (6, 3, 4)]:
            replacement = {n:size, k:bound, h:stride}
            concrete = view.xreplace(replacement)
            expected_raw = slice(*(None if item is None else int(sx.sympify(item).xreplace(replacement)) for item in (raw.start, raw.stop, raw.step)))
            expected = list(range(size))[expected_raw]
            assert concrete.shape == (len(expected), 4)
            assert [concrete[i, 0] for i in range(len(expected))] == [sx.MatrixSymbol('A', size, 4)[i, 0] for i in expected]
            assert concrete.func(*concrete.args) == concrete
            assert sx.sympify(sx.srepr(concrete)) == concrete
    with pytest.raises(NotImplementedError, match='sign'):
        A[::k, :]


def test_equivalent_normalized_ranges_have_identical_hashes_and_args():
    A = sx.MatrixSymbol('A', 5, 6)
    for left, right in [(A[:5:2, :], A[:6:2, :]), (A[100:200, :], A[2:2:-1, :])]:
        assert left == right
        assert hash(left) == hash(right)
        assert left.args == right.args
        assert len({left, right}) == 1


def test_symbolic_scalar_index_keeps_unit_axis_and_checks_substitution():
    n = sx.Symbol('n', integer=True, positive=True)
    i = sx.Symbol('i', integer=True)
    A = sx.MatrixSymbol('A', n, 4)
    view = A[i, ::-1]
    assert view.shape == (1, 4)
    concrete = view.xreplace({n:5, i:-1})
    assert concrete[0, 0] == sx.MatrixSymbol('A', 5, 4)[4, 3]
    with pytest.raises(IndexError):
        view.xreplace({n:5, i:5})


def test_arbitrary_precision_bounds_steps_and_scalar_out_of_range():
    parent = numeric_parent(5, 6)
    huge = 10**100
    for row in [slice(-huge, huge, 2), slice(huge, -huge, -2), slice(None, None, huge), slice(None, None, -huge)]:
        assert_numeric_view(parent[row, :], (5, 6), row, slice(None))
    for index in [huge, -huge]:
        with pytest.raises(IndexError):
            parent[index, :]
        with pytest.raises(IndexError):
            parent[:, index]


@pytest.mark.parametrize('protocol', range(pickle.HIGHEST_PROTOCOL+1))
def test_public_normalized_slice_tuple_reconstructs_after_standalone_pickle(protocol):
    A = sx.MatrixSymbol('A', 5, 6)
    view = A[::-2, 1::2]
    rows = pickle.loads(pickle.dumps(view.rowslice, protocol=protocol))
    cols = pickle.loads(pickle.dumps(view.colslice, protocol=protocol))
    assert isinstance(rows, tuple) and rows == view.rowslice
    assert rows[1] is None
    assert sx.MatrixSlice(A, rows, cols) == view
