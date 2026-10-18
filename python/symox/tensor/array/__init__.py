"""Concrete multidimensional arrays and tensor operations on native scalars."""
from ._array import (NDimArray, DenseNDimArray, SparseNDimArray,
    MutableNDimArray, ImmutableNDimArray, MutableDenseNDimArray,
    ImmutableDenseNDimArray, MutableSparseNDimArray, ImmutableSparseNDimArray, Array)
from .arrayop import tensorproduct, tensorcontraction, tensordiagonal, permutedims, derive_by_array
from .array_comprehension import ArrayComprehension,ArrayComprehensionMap
from .kind import ArrayKind

__all__=['Array','NDimArray','DenseNDimArray','SparseNDimArray','MutableNDimArray',
    'ImmutableNDimArray','MutableDenseNDimArray','ImmutableDenseNDimArray',
    'MutableSparseNDimArray','ImmutableSparseNDimArray','tensorproduct',
    'tensorcontraction','tensordiagonal','permutedims','derive_by_array']
__all__ += ['ArrayComprehension','ArrayComprehensionMap','ArrayKind']
