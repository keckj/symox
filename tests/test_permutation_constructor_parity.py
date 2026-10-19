"""Cycle notation must survive upstream srepr and overlapping composition."""
import itertools
import random

import pytest
from sympy.combinatorics import Cycle as ReferenceCycle
from sympy.combinatorics import Permutation as ReferencePermutation
from symox.combinatorics import Cycle, Permutation


@pytest.mark.parametrize('args', [(), (-2,), (-1,), (0,), (True,), (4,),
    (0, 1), (3, 1, 2), ([],), ([1, 0],), ([[0, 1], [1, 2]],),
    ([[3], [], [0, 2]],)])
@pytest.mark.parametrize('size', [None, 0, 1, 5])
def test_constructor_forms(args, size):
    try:
        expected = ReferencePermutation(*args, size=size)
    except (ValueError, TypeError) as error:
        with pytest.raises(type(error)):
            Permutation(*args, size=size)
    else:
        assert Permutation(*args, size=size).array_form == expected.array_form


def test_overlapping_cycles_and_chained_calls():
    rng = random.Random(1729)
    for degree in range(1, 20):
        cycles = [rng.sample(range(degree), rng.randrange(degree + 1))
                  for _ in range(12)]
        actual = Permutation(cycles)
        expected = ReferencePermutation(cycles)
        assert actual.array_form == expected.array_form
        assert Permutation(Cycle(*cycles[0]), size=degree).array_form == ReferencePermutation(
            ReferenceCycle(*cycles[0]), size=degree).array_form
        for cycle in cycles:
            if len(cycle) == 1:
                assert actual(cycle[0]) == expected(cycle[0])
            else:
                actual, expected = actual(*cycle), expected(*cycle)
                assert actual.array_form == expected.array_form


def test_application_to_indices_and_sequences():
    for images in itertools.permutations(range(4)):
        actual, expected = Permutation(images), ReferencePermutation(images)
        for sequence in ['abcd', [10, 20, 30, 40], (None, False, 0, 1)]:
            assert actual(sequence) == expected(sequence)
        for invalid in [-2, -1, 5, 0.5, [1], 'abc']:
            with pytest.raises(TypeError):
                actual(invalid)


@pytest.mark.parametrize('size', [None, -1, 0, 2, 6])
def test_cycle_list_preserves_degree(size):
    assert Cycle(2, 3).list(size) == ReferenceCycle(2, 3).list(size)
    assert Cycle(2, 3)(8).list(size) == ReferenceCycle(2, 3)(8).list(size)


def test_empty_cycle_requires_explicit_size():
    with pytest.raises(ValueError):
        Cycle().list()
    with pytest.raises(ValueError):
        Permutation(Cycle())
