"""LaTeX printing and lambdify of N-dimensional arrays, compared with SymPy."""
import numpy as np
import pytest
import sympy as sp
import symox as sx

_ARRAYS = {
    'matrix': lambda m, x, y: m.Array([[x, y], [0, x*y]]),
    'vector': lambda m, x, y: m.Array([x, y**2]),
    'scalar': lambda m, x, y: m.Array(x),
    'rank three': lambda m, x, y: m.Array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]]),
    'wide': lambda m, x, y: m.Array(list(range(12))),
    'mutable': lambda m, x, y: m.MutableDenseNDimArray([[x, 1/y]]),
    'sparse': lambda m, x, y: m.ImmutableSparseNDimArray([x, 0, y]),
    'subclass': lambda m, x, y: type('Sub', (m.ImmutableDenseNDimArray,), {})([[x], [y]]),
    'one by one': lambda m, x, y: m.Array([[x]]),
    'product': lambda m, x, y: x*m.Array([x, y]),
}
_OPTIONS = [{}, {'mode': 'inline'}, {'mat_delim': '('}, {'mat_delim': ''}, {'mat_str': 'array'}]


@pytest.mark.parametrize('options', _OPTIONS)
@pytest.mark.parametrize('name', list(_ARRAYS))
def test_array_latex_matches_sympy(name, options):
    def run(lib):
        x, y = lib.symbols('x y')
        return lib.latex(_ARRAYS[name](lib, x, y), **options)
    assert run(sx) == run(sp)


def test_empty_array_latex_raises_like_sympy():
    for lib in (sx, sp):
        with pytest.raises(IndexError):
            lib.latex(lib.Array([], (0,)))


@pytest.mark.parametrize('modules', [None, np, 'numpy', ['numpy'], 'scipy'])
@pytest.mark.parametrize('build, args', [
    (lambda m, x, y: m.Array([[x, y], [1, x*y]]), (1.5, 2)),
    (lambda m, x, y: m.Array([m.sin(x), y]), (1, 2)),
    (lambda m, x, y: m.ImmutableSparseNDimArray([x, 0, y]), (1, 2)),
    (lambda m, x, y: m.Array(x*y), (3, 2)),
    (lambda m, x, y: m.Array([], (2, 0)), (3, 2)),
    (lambda m, x, y: m.Array([x, x*y]), (np.array([1., 2.]), np.array([3., 4.]))),
])
def test_array_lambdify_values_match_sympy(modules, build, args):
    def run(lib):
        x, y = lib.symbols('x y')
        return lib.lambdify([x, y], build(lib, x, y), modules)(*args)
    actual, expected = run(sx), run(sp)
    assert isinstance(actual, np.ndarray) and actual.shape == expected.shape
    np.testing.assert_allclose(actual, expected)


@pytest.mark.parametrize('modules', ['math', 'mpmath'])
def test_array_lambdify_without_numpy_is_rejected(modules):
    x, y = sx.symbols('x y')
    with pytest.raises(NotImplementedError):
        sx.lambdify([x, y], sx.Array([[x, y], [1, x*y]]), modules)(1, 2)


def test_mutable_array_lambdify_is_rejected_like_sympy():
    for lib in (sx, sp):
        x, y = lib.symbols('x y')
        with pytest.raises(NotImplementedError):
            lib.lambdify([x, y], lib.MutableDenseNDimArray([x, y]))(3, 2)
