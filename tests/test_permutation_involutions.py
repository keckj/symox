"""Exact involution reduction covers every signed integer exponent."""
import pytest
from symox.combinatorics import Permutation
from sympy.combinatorics import Permutation as Reference


@pytest.mark.parametrize('size', [0,1,2,3,20,101])
@pytest.mark.parametrize('exponent', [-2**100-1,-17,-8,0,8,9,17,2**100,2**100+1])
def test_disjoint_transpositions_and_fixed_points_match_reference(size,exponent):
    images = list(range(size))
    for i in range(0,size-1,3):
        images[i],images[i+1] = images[i+1],images[i]
    value = Permutation(images)
    value.rank();value.parity();value.cycle_structure
    result = value**exponent
    assert result.array_form == (Reference(images)**exponent).array_form
    assert result.size == size
    assert result.rank() == Reference(result.array_form).rank()


def test_fixed_prefix_does_not_hide_a_longer_cycle():
    images = list(range(12))+[13,14,15,12,17,16]
    for exponent in [-17,8,9,17,2**100+1]:
        assert (Permutation(images)**exponent).array_form == (Reference(images)**exponent).array_form
