"""Group enumeration, conjugation and uniform sampling used by Haarpy."""
import inspect
import itertools
import math

import pytest
from symox.combinatorics import Permutation, PermutationGroup
from symox.combinatorics.named_groups import SymmetricGroup, AlternatingGroup, DihedralGroup, CyclicGroup
from sympy.combinatorics import named_groups as reference


@pytest.mark.parametrize('factory,size', [(SymmetricGroup, 1), (SymmetricGroup, 4),
    (AlternatingGroup, 5), (DihedralGroup, 5), (CyclicGroup, 7)])
def test_generators_and_conjugation_match_reference(factory, size):
    native, expected = factory(size), getattr(reference, factory.__name__)(size)
    normalize = lambda values: {tuple(value.array_form) for value in values}
    for method in ['coset', 'dimino']:
        actual = list(native.generate(method=method))
        assert len(actual) == expected.order()
        assert normalize(actual) == normalize(expected.generate(method=method))
        assert {tuple(v) for v in native.generate(method=method, af=True)} == normalize(actual)
    # Dimino order depends on the chosen generators; named groups can use
    # different (equivalent) generating sets.
    from sympy.combinatorics import Permutation as RP, PermutationGroup as RG
    same_generators = RG([RP(g.array_form) for g in native.generators])
    assert list(native.generate_dimino(af=True)) == list(same_generators.generate_dimino(af=True))
    classes = native.conjugacy_classes()
    assert {frozenset(normalize(c)) for c in classes} == {
        frozenset(normalize(c)) for c in expected.conjugacy_classes()}
    assert classes[0] == {native.identity}
    assert sum(map(len, classes)) == native.order()
    for a, b in zip(native.generators, same_generators.generators):
        assert normalize(native.conjugacy_class(a)) == normalize(same_generators.conjugacy_class(b))
    for method in ['generate', 'generate_dimino', 'random', 'conjugacy_class', 'conjugacy_classes']:
        assert inspect.signature(getattr(native, method)) == inspect.signature(getattr(expected, method))
    with pytest.raises(NotImplementedError): native.generate('unknown')


def test_conjugacy_class_accepts_elements_outside_group():
    group = CyclicGroup(4)
    element = Permutation(0, 1, size=4)
    expected = reference.CyclicGroup(4)
    from sympy.combinatorics import Permutation as RP
    assert {tuple(v.array_form) for v in group.conjugacy_class(element)} == {
        tuple(v.array_form) for v in expected.conjugacy_class(RP(0, 1, size=4))}


def test_transversal_sampling_is_uniform_without_enumerating_group(monkeypatch):
    import symox.combinatorics.perm_groups as module
    # As in the reference, random() unranks a uniform rank through the
    # basic transversals, so every rank gives a distinct group element.
    group = SymmetricGroup(4)
    seen = set()
    for rank in range(group.order()):
        monkeypatch.setattr(module, 'randrange', lambda n, rank=rank: rank)
        value = group.random()
        assert group.contains(value)
        seen.add(tuple(value.array_form))
    assert len(seen) == group.order() == math.prod(map(len, group.basic_transversals))
    large = SymmetricGroup(30)
    monkeypatch.setattr(module, 'randrange', lambda n: n - 1)
    def forbidden(*args, **kwargs):
        raise AssertionError('Sampling must not enumerate the factorial-size group')
    monkeypatch.setattr(large, 'generate_schreier_sims', forbidden)
    monkeypatch.setattr(large, '_enumerate', forbidden)
    assert large.contains(large.random())
    assert large.contains(Permutation(large.random(af=True)))
