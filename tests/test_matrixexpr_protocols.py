"""Preserve matrix types, shapes and factor order through Basic protocols."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def pair():
    n = sx.Symbol('n', positive=True, integer=True)
    x, y = sx.symbols('x y')
    a, b = sx.MatrixSymbol('A', n, n), sx.MatrixSymbol('B', n, n)
    return x, y, a, b


def convert(value):
    if hasattr(value, 'to_sympy'):
        return value.to_sympy()
    if isinstance(value, dict):
        return {convert(k): convert(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return type(value)(convert(v) for v in value)
    return value


@pytest.mark.parametrize('index', range(9))
@pytest.mark.parametrize('method', ['atoms', 'as_coeff_add', 'as_numer_denom',
                                   'as_ordered_factors', 'as_powers_dict', 'args_cnc',
                                   'copy', 'getO', 'getn', 'removeO', 'primitive'])
def test_structural_views(index, method):
    x, y, a, b = pair()
    e = [a, a.T, a.I, a**2, a*b, 2*x*a, 2*a+3*b,
         sx.Identity(a.rows), sx.ZeroMatrix(a.rows, a.cols)][index]
    assert convert(getattr(e, method)()) == getattr(e.to_sympy(), method)()


@pytest.mark.parametrize('index', range(5))
def test_traversal_and_free_queries(index):
    x, y, a, b = pair()
    e = [a, a*b, 2*x*a+b, a.I, sx.Identity(2)][index]
    reference = e.to_sympy()
    assert convert(e.atoms(sx.Symbol)) == reference.atoms(sp.Symbol)
    assert convert(e.atoms(sx.MatrixSymbol)) == reference.atoms(sp.MatrixSymbol)
    assert convert(e.find(sx.MatrixSymbol)) == reference.find(sp.MatrixSymbol)
    assert convert(e.find(sx.MatrixSymbol, group=True)) == reference.find(sp.MatrixSymbol, group=True)
    assert e.count(sx.MatrixSymbol) == reference.count(sp.MatrixSymbol)
    for needle in [x, y, a, b, a.rows, sx.Integer(2)]:
        assert e.has(needle) == reference.has(convert(needle))
        assert e.has_free(needle) == reference.has_free(convert(needle))
        assert e.has_xfree({needle}) == reference.has_xfree({convert(needle)})
    assert e.has(sx.MatrixSymbol) == reference.has(sp.MatrixSymbol)
    with pytest.raises(TypeError):
        e.has_xfree([x])


@pytest.mark.parametrize('index', range(7))
@pytest.mark.parametrize('method', ['evalf', 'n', 'simplify', 'trigsimp', 'cancel',
                                   'factor', 'together', 'powsimp', 'radsimp'])
def test_scalar_transformations(index, method):
    x, y, a, b = pair()
    e = [a, a.I, a**2, sx.Identity(2),
         sx.sin(x)**2*a+sx.cos(x)**2*a,
         ((x*x-1)/(x-1))*a*b,
         (sx.Rational(1, 3)*x)*a+b][index]
    actual = getattr(e, method)()
    # Several scalar simplifiers in SymPy 1.14 raise on MatrixExpr. In
    # those cases compare each explicit entry with the scalar operation.
    try:
        if method == 'radsimp':
            raise TypeError('SymPy radsimp can return a scalar Mul containing matrix operands')
        expected = getattr(e.to_sympy(), method)()
    except (TypeError, AttributeError):
        substitutions = {convert(a.rows): 2}
        expected = e.to_sympy().subs(substitutions).as_explicit().applyfunc(
            lambda entry: getattr(entry, method)())
        actual_explicit = actual.to_sympy().subs(substitutions).as_explicit()
        assert (actual_explicit-expected).applyfunc(sp.simplify) == sp.zeros(2)
        return
    assert isinstance(actual, sx.MatrixExpr)
    assert actual.shape == e.shape
    assert sp.simplify(actual.to_sympy()-expected) == sp.ZeroMatrix(*expected.shape)


@pytest.mark.parametrize('index', range(5))
def test_real_imag_and_positive_content(index):
    x, y, a, b = pair()
    e = [a, a.T, 2*x*a, sx.Rational(2, 3)*a*b, 2*a+4*b][index]
    re, im = e.as_real_imag()
    reference_re, reference_im = e.to_sympy().as_real_imag()
    assert sp.simplify(convert(re)-reference_re) == sp.ZeroMatrix(*e.shape)
    assert sp.simplify(convert(im)-reference_im) == sp.ZeroMatrix(*e.shape)
    content, primitive = e.as_content_primitive()
    assert content.is_positive
    assert (content*primitive).doit() == e.doit()


def test_structural_matching_and_repeated_wildcards():
    x, y, a, b = pair()
    w = sx.Wild('w')
    pattern = sx.MatMul(w, a, b)
    assert pattern.matches(3*x*a*b) == {w: 3*x}
    assert (3*x*a*b).match(pattern) == {w: 3*x}
    assert pattern.matches(3*x*b*a) is None
    assert a.matches(a) == {}
    assert a.matches(b) is None
    pattern = sx.MatAdd(w*a, w*b)
    assert pattern.matches(3*a+3*b) == {w: sx.Integer(3)}
    assert pattern.matches(3*a+4*b) is None


def test_dummy_equivalence_and_scoped_scalar_queries():
    x, y, a, b = pair()
    i, j, p, q = [sx.Dummy(v, integer=True) for v in 'ijpq']
    left = sx.FunctionMatrix(2, 2, sx.Lambda((i, j), i+2*j+x))
    right = sx.FunctionMatrix(2, 2, sx.Lambda((p, q), p+2*q+x))
    assert left.dummy_eq(right)
    assert not left.dummy_eq(sx.FunctionMatrix(2, 2, sx.Lambda((p, q), p+2*q+y)))
    assert not left.has_free(i) and left.has_free(x)
    assert not left.has_xfree({i}) and left.has_xfree({x})
    assert left.as_dummy().as_explicit() == left.as_explicit()


@pytest.mark.parametrize('i,j', [(-3, 0), (-2, 0), (-1, 0), (0, 0), (1, 2),
                               (2, 0), (0, 3), (sx.Symbol('i'), 0)])
def test_index_validation(i, j):
    a = sx.MatrixSymbol('A', 2, 3)
    assert a.valid_index(i, j) == a.to_sympy().valid_index(convert(i), convert(j))


def test_structural_equality_and_comparison():
    _, _, a, b = pair()
    for left in [a, b, a.T, a*b, b*a]:
        assert left.is_same(left)
        assert not left.is_same(left.to_sympy())
        for right in [a, b, a.T, a*b, b*a]:
            assert left.compare(right) == left.to_sympy().compare(right.to_sympy())


def test_native_matrix_protocols_without_sympy():
    code = '''
import sys
sys.modules['sympy']=None
import symox as s
x=s.Symbol('x'); a=s.MatrixSymbol('A',2,2); b=s.MatrixSymbol('B',2,2)
e=(s.sin(x)**2+s.cos(x)**2)*a*b
assert e.simplify() == a*b
assert e.find(s.MatrixSymbol) == {a,b}
assert (2*x*a).has_free(x)
assert a.I.evalf().shape == (2,2)
assert (2*a+4*b).as_content_primitive()[0] == 2
'''
    subprocess.run([sys.executable, '-c', code], check=True, timeout=30)
