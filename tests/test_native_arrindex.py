"""Native (Rust) array_derive, convert_indexed_to_array and
convert_array_to_indexed (docs/algorithmic-parity.md, section arrindex).

The differential corpus of audits/performance/2026-10-30-arrindex is evaluated
with symox and compared with the SymPy 1.14 results recorded in
sympy_expected.json. The cases listed in CORE_ORDER differ from SymPy only by
the order of the factors of a scalar Mul/Add (symox core canonical ordering,
a preexisting gap outside these conversions); their symox result is pinned.
"""
import importlib.util
import json
import pathlib

import pytest

import symox as S
from symox import IndexedBase, KroneckerDelta, MatrixSymbol, Sum, symbols
from symox.tensor.array.expressions import (ArrayAdd, ArrayContraction, ArrayDiagonal, ArrayElement, ArraySymbol,
                                            ArrayTensorProduct, OneArray)
from symox.tensor.array.expressions.arrayexpr_derivatives import array_derive
from symox.tensor.array.expressions.from_array_to_indexed import convert_array_to_indexed
from symox.tensor.array.expressions.from_indexed_to_array import convert_indexed_to_array

AUDIT = pathlib.Path(__file__).resolve().parent.parent / 'audits' / 'performance' / '2026-10-30-arrindex'


def _corpus():
    spec = importlib.util.spec_from_file_location('_arrindex_corpus', AUDIT / 'corpus.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CORPUS = _corpus()
EXPECTED = json.loads((AUDIT / 'sympy_expected.json').read_text())

CORE_ORDER = {
    ('indexed', 'add_bcast'): 'ArrayElement(ArrayAdd(ArrayTensorProduct(MatrixSymbol:v, OneArray(d,)), PermuteDims('
                              'ArrayTensorProduct(MatrixSymbol:M, OneArray(1,)), [0, 2, 1])), [i, 0, j])',
    ('indexed', 'int_mul'): 'ArrayElement(ArrayTensorProduct(MatrixSymbol:N, MatrixSymbol:M), [0, j, i, 0])',
    ('indexed', 'int_mul_first'): 'ArrayElement(ArrayTensorProduct(MatrixSymbol:N, MatrixSymbol:M), [0, j, i, 0])',
    ('indexed', 'u_int_mul'): '(ArrayDiagonal(ArrayTensorProduct(IndexedBase:IA, IndexedBase:IB), (0, 2)), '
                              '(Zero:0, Symbol:i))',
    ('to_indexed', 'numbers'): 'Mul:C[2]*B[0, 1]',
    ('to_indexed', 'tp'): 'Mul:C[l]*B[i, j]',
}

CASES = [(s, n) for s, cases in CORPUS.SECTIONS.items() for n in sorted(cases)]


@pytest.mark.parametrize('section,name', CASES, ids=['%s-%s' % c for c in CASES])
def test_corpus_matches_sympy(section, name):
    got = CORPUS.evaluate('symox', section, name)
    assert got == CORE_ORDER.get((section, name), EXPECTED[section][name])


def test_native_entry_points_are_used(monkeypatch):
    from symox import _core
    calls = []
    for name in ('_arrayexpr_derive', '_arrayexpr_from_indexed', '_arrayexpr_to_indexed'):
        original = getattr(_core, name)

        def wrapper(*args, _original=original, _name=name):
            calls.append(_name)
            return _original(*args)
        monkeypatch.setattr(_core, name, wrapper)
    k, i, j = symbols('k i j')
    X = MatrixSymbol('X', k, k)
    array_derive(X, X)
    A = ArraySymbol('A', (k, k))
    convert_indexed_to_array(Sum(A[i, j], (j, 0, k - 1)))
    convert_array_to_indexed(ArrayContraction(A, (0, 1)), [])
    assert calls == ['_arrayexpr_derive', '_arrayexpr_from_indexed', '_arrayexpr_to_indexed']


def test_same_named_arrays_keep_their_parent():
    i, j, k, l = symbols('i j k l')
    small = ArraySymbol('A', (k,))
    big = ArraySymbol('A', (3, 3))
    B = ArraySymbol('B', (3, 3))
    assert convert_indexed_to_array(small[i]) == small
    O2 = OneArray(3, 3)
    assert convert_indexed_to_array(big[i, j] + B[k, l]) == ArrayAdd(ArrayTensorProduct(big, O2),
                                                                     ArrayTensorProduct(O2, B))
    # Elements of distinct arrays are distinct (SymPy compares the arrays).
    assert ArrayElement(small, (i,))._symox_() != ArrayElement(ArraySymbol('A', (3,)), (i,))._symox_()
    assert ArrayElement(big, (i, j))._symox_() == ArrayElement(ArraySymbol('A', (3, 3)), (i, j))._symox_()
    assert S.diff(big[i, j], big[k, l]) == KroneckerDelta(i, k)*KroneckerDelta(j, l)


def test_symbolic_dimension_index():
    # SymPy's MatrixExpr.valid_index leaves `k < k` unevaluated for a plain
    # symbol, so `M[i, k]` is a valid entry of a k x k matrix.
    i, j, k, l = symbols('i j k l')
    M = MatrixSymbol('M', k, k)
    N = MatrixSymbol('N', k, k)
    assert str(M[i, k]) == 'M[i, k]'
    assert convert_indexed_to_array(M[i, j]*N[k, l], [i, j, k, l]) == ArrayTensorProduct(M, N)
    n = symbols('n', integer=True, positive=True)
    with pytest.raises(IndexError):
        MatrixSymbol('Z', n, n)[n, 0]


def test_errors():
    i, j, k = symbols('i j k')
    C = ArraySymbol('C', (3, 4))
    D = ArraySymbol('D', (4, 3))
    with pytest.raises(ValueError):
        convert_indexed_to_array(C[i, j] + D[i, j])
    IA = IndexedBase('IA')
    # SymPy raises TypeError on unshaped bases in these two cases.
    with pytest.raises(TypeError):
        convert_indexed_to_array(IA[i, j] + IA[j, k])
    with pytest.raises(TypeError):
        convert_indexed_to_array(Sum(IA[i, j]*IndexedBase('IB')[j], (j, 0, 3)), [i])
    with pytest.raises(NotImplementedError):
        array_derive(object(), ArraySymbol("A", (3,)))


def test_unshaped_indexed_base():
    i, j, k = symbols('i j k')
    IA, IB = IndexedBase('IA'), IndexedBase('IB')
    result = convert_indexed_to_array(IA[i, j]*IB[j, k])
    assert result == ArrayDiagonal(ArrayTensorProduct(IA, IB), (1, 2))
    assert result.shape is None
    assert convert_indexed_to_array(Sum(IA[i, j]*IB[j, k], (j, 0, 3))).shape is None
