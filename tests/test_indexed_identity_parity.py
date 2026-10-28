"""Indexed scalar atoms preserve bounds and shape across shared printer state."""
import pickle

import pytest

import symox as s
from symox.tensor.indexed import IndexException


def test_unbounded_reference_does_not_replace_later_bounded_reference():
    base = s.IndexedBase('identity_y')
    unbounded = base[s.Idx('identity_i')]._symox_()
    bounded = base[s.Idx('identity_i', 3)]._symox_()
    assert unbounded != bounded
    assert s.jscode(bounded + 1) == 'identity_y[identity_i] + 1'
    with pytest.raises(IndexException):
        s.jscode(unbounded + 1)
    assert s.jscode(bounded + 1) == 'identity_y[identity_i] + 1'


@pytest.mark.parametrize('explicit_shape', [False, True])
def test_same_names_with_two_shapes_coexist(explicit_shape):
    def indexed(size):
        base = s.IndexedBase('identity_A', shape=(size,) if explicit_shape else None)
        return base[s.Idx('identity_j', size)]._symox_()
    first, second = indexed(3), indexed(5)
    assert first != second
    assert len({first, second}) == 2
    assert (first + second).subs(first, 0) == second
    from symox.codegen.ast import _indexed_reference
    assert tuple(_indexed_reference(first).shape) == (s.Integer(3),)
    assert tuple(_indexed_reference(second).shape) == (s.Integer(5),)
    for value, size in [(first, 3), (second, 5), (first, 3)]:
        code = s.ccode(value + 1, assign_to=s.IndexedBase('identity_out')[s.Idx('identity_j', size)], contract=True)
        assert 'identity_j<%s' % size in code.replace(' ', '')


def test_same_metadata_from_separate_allocations_compares_equal():
    def indexed():
        return s.IndexedBase('identity_equal', shape=(7,))[s.Idx('identity_equal_i', (2, 8))]._symox_()
    first, second = indexed(), indexed()
    assert first == second
    assert hash(first) == hash(second)
    assert len({first, second}) == 1
    assert first + second == 2 * first
    assert first != s.Symbol(str(first))


def test_index_bounds_are_part_of_identity_even_with_same_dimension():
    base = s.IndexedBase('identity_shift')
    first = base[s.Idx('identity_shift_i', (0, 2))]._symox_()
    second = base[s.Idx('identity_shift_i', (1, 3))]._symox_()
    assert first != second
    from symox.codegen.ast import _indexed_reference
    assert _indexed_reference(first).indices[0].lower == 0
    assert _indexed_reference(second).indices[0].lower == 1


@pytest.mark.parametrize('compound', [False, True])
def test_pickle_retains_metadata_without_process_registry(monkeypatch, compound):
    import symox.codegen.ast as ast
    import symox._public_matrix as public
    first = s.IndexedBase('identity_pickle')[s.Idx('identity_pickle_i', 4)]._symox_()
    expression = first + 1 if compound else first
    encoded = pickle.dumps(expression)
    monkeypatch.setattr(ast, '_INDEXED_SYMBOLS', {})
    monkeypatch.setattr(public, '_indexed', {})
    decoded = pickle.loads(encoded)
    assert decoded == expression
    assert hash(decoded) == hash(expression)
    assert s.jscode(decoded) == 'identity_pickle[identity_pickle_i]' + (' + 1' if compound else '')
    assert tuple(ast._indexed_reference(first).shape) == (s.Integer(4),)


def test_explicit_sympy_roundtrip_preserves_indexed_metadata():
    sp = pytest.importorskip('sympy')
    index = s.Idx('identity_roundtrip_i', (2, 6))
    value = s.IndexedBase('identity_roundtrip', shape=(7,))[index]._symox_()
    foreign = s.to_sympy(value + 1)
    expected = sp.IndexedBase('identity_roundtrip', shape=(7,))[sp.Idx('identity_roundtrip_i', (2, 6))] + 1
    assert foreign == expected
    assert s.from_sympy(foreign) == value + 1


def test_unshaped_multidimensional_metadata_preserved():
    from symox.codegen.ast import _indexed_reference
    value = s.IndexedBase('identity_unshaped')[s.Idx('identity_u'), s.Idx('identity_v')]._symox_()
    restored = _indexed_reference(pickle.loads(pickle.dumps(value)))
    assert restored.base.shape is None
    assert all(i.lower is None for i in restored.indices)
