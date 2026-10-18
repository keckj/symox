"""Concrete multidimensional arrays and tensor operations on native scalars."""
from ._array import (NDimArray, DenseNDimArray, SparseNDimArray,
    MutableNDimArray, ImmutableNDimArray, MutableDenseNDimArray,
    ImmutableDenseNDimArray, MutableSparseNDimArray, ImmutableSparseNDimArray, Array)
from .arrayop import tensorproduct, tensorcontraction, tensordiagonal, permutedims, derive_by_array

__all__=['Array','NDimArray','DenseNDimArray','SparseNDimArray','MutableNDimArray',
    'ImmutableNDimArray','MutableDenseNDimArray','ImmutableDenseNDimArray',
    'MutableSparseNDimArray','ImmutableSparseNDimArray','tensorproduct',
    'tensorcontraction','tensordiagonal','permutedims','derive_by_array']
