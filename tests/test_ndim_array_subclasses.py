"""User subclasses of the N-dimensional array classes, compared with SymPy."""
import pytest
import sympy as sp
import symox as sx

_BASES = ['ImmutableDenseNDimArray', 'MutableDenseNDimArray',
          'ImmutableSparseNDimArray', 'MutableSparseNDimArray']

_OPERATIONS = {
    'construct': lambda a, m, x, y: a,
    'add': lambda a, m, x, y: a + a,
    'sub': lambda a, m, x, y: a - a,
    'mul': lambda a, m, x, y: a * 2,
    'rmul': lambda a, m, x, y: 2 * a,
    'mul zero': lambda a, m, x, y: a * 0,
    'div': lambda a, m, x, y: a / 2,
    'neg': lambda a, m, x, y: -a,
    'slice': lambda a, m, x, y: a[0, :],
    'row': lambda a, m, x, y: a[0],
    'item': lambda a, m, x, y: a[1, 1],
    'flat item': lambda a, m, x, y: a[3],
    'negative index': lambda a, m, x, y: a[-1, -1],
    'applyfunc': lambda a, m, x, y: a.applyfunc(lambda e: e + 1),
    'reshape': lambda a, m, x, y: a.reshape(4),
    'diff': lambda a, m, x, y: a.diff(x),
    'diff function': lambda a, m, x, y: m.diff(a, x),
    'transpose': lambda a, m, x, y: a.transpose(),
    'tolist': lambda a, m, x, y: a.tolist(),
    'tensorcontraction': lambda a, m, x, y: m.tensorcontraction(a, (0, 1)),
    'derive_by_array': lambda a, m, x, y: m.derive_by_array(a, [x, y]),
    'permutedims': lambda a, m, x, y: m.permutedims(a, (1, 0)),
    'tensorproduct type': lambda a, m, x, y: type(m.tensorproduct(a, [1, 2])).__name__,
    'equal': lambda a, m, x, y: a == type(a)([[x, y], [0, x*y]]),
    'equal list': lambda a, m, x, y: a == [[x, y], [0, x*y]],
    'len': lambda a, m, x, y: len(a),
    'iter': lambda a, m, x, y: list(a),
    'zeros': lambda a, m, x, y: type(a).zeros(2, 2),
    'rank': lambda a, m, x, y: a.rank(),
    'str': lambda a, m, x, y: str(a),
    'repr': lambda a, m, x, y: repr(a),
    'latex': lambda a, m, x, y: m.latex(a),
    'sympify type': lambda a, m, x, y: type(m.sympify(a)).__name__,
    'nested arrays': lambda a, m, x, y: type(a)([type(a)([1, 2]), type(a)([3, 4])]),
    'flat with shape': lambda a, m, x, y: type(a)([1, 2, 3, 4], (2, 2)),
    'from matrix': lambda a, m, x, y: type(a)(m.Matrix([[1, x], [2, 3]])),
    'multiply arrays': lambda a, m, x, y: a * a,
    'add list': lambda a, m, x, y: a + [[1, 2], [3, 4]],
    'out of bounds': lambda a, m, x, y: a[2, 0],
    'too many indices': lambda a, m, x, y: a[0, 0, 0],
    'ragged': lambda a, m, x, y: type(a)([[1, 2], [3]]),
}


def _run(lib, base, operation):
    cls = type('Sub', (getattr(lib, base),), {})
    x, y = lib.symbols('x y')
    array = cls([[x, y], [0, x*y]])
    try:
        result = _OPERATIONS[operation](array, lib, x, y)
    except Exception as error:
        return 'raises ' + type(error).__name__
    name = type(result).__name__
    if name in {'MutableDenseMatrix', 'Matrix'}:
        name = 'Matrix'
    return name, str(result)


@pytest.mark.parametrize('operation', list(_OPERATIONS))
@pytest.mark.parametrize('base', _BASES)
def test_array_subclass_results_match_sympy(base, operation):
    assert _run(sx, base, operation) == _run(sp, base, operation)


def test_subclass_with_extra_constructor_arguments():
    # SymPy calls ``__init__`` after ``__new__``; derived arrays are rebuilt
    # through the public constructor, which such a subclass rejects.
    for lib in (sx, sp):
        x = lib.Symbol('x')

        class Tagged(lib.ImmutableDenseNDimArray):
            def __new__(cls, iterable, a, b, shape=None, **kwargs):
                return super().__new__(cls, iterable, shape, **kwargs)

            def __init__(self, iterable, a, b, shape=None, **kwargs):
                self.a = a

        tagged = Tagged([[x, 1], [2, 3]], 1, 2)
        assert type(tagged) is Tagged and tagged.a == 1
        assert list(tagged._array) == [x, 1, 2, 3]
        assert tagged.free_symbols == {x}
        assert tagged.tomatrix() == lib.Matrix([[x, 1], [2, 3]])
        for operation in (lambda: tagged + tagged, lambda: tagged * 2, lambda: -tagged,
                          lambda: tagged.diff(x), lambda: tagged.applyfunc(lambda e: e),
                          lambda: tagged.transpose(), lambda: tagged.reshape(4),
                          lambda: tagged.subs(x, 1), lambda: tagged[0, :]):
            with pytest.raises(TypeError):
                operation()


@pytest.mark.parametrize('call', [
    lambda C, x: C._handle_ndarray_creation_inputs([[1, 2], [3, 4]]),
    lambda C, x: C._handle_ndarray_creation_inputs([1, 2, 3, 4], (2, 2)),
    lambda C, x: C._handle_ndarray_creation_inputs(x),
    lambda C, x: C._scan_iterable_shape([[1, 2], [3, 4]]),
    lambda C, x: C._new([1, 2], (2,)),
    lambda C, x: C([[1, 2], [3, 4]])._get_tuple_index(3),
    lambda C, x: C([[1, 2], [3, 4]])._check_index_for_getitem(1),
    lambda C, x: C([[x, 2], [3, 4]])._loop_size,
    lambda C, x: C([[x, 2], [3, 4]])._rank,
])
def test_array_protocol_helpers_match_sympy(call):
    def run(lib):
        result = call(lib.ImmutableDenseNDimArray, lib.Symbol('x'))
        return str(result)
    assert run(sx) == run(sp)


@pytest.mark.parametrize('kind', ['ImmutableSparseNDimArray', 'MutableSparseNDimArray'])
def test_dense_conversion_of_sparse_arrays_keeps_the_values(kind):
    # SymPy 1.14 builds a dense array from the stored flat positions of a
    # sparse one (ImmutableDenseNDimArray(ImmutableSparseNDimArray([x, y]))
    # is [0, 1]); symox keeps the entries, which is the documented meaning.
    x, y = sx.symbols('x y')
    sparse = getattr(sx, kind)([[x, y], [0, x*y]])
    dense = sx.ImmutableDenseNDimArray(sparse)
    assert dense.shape == (2, 2) and dense.tolist() == [[x, y], [0, x*y]]
    assert sx.derive_by_array(x*y, getattr(sx, kind)([x, y])) == sx.Array([y, x])
