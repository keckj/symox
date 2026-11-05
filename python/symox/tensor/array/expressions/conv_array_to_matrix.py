"""Deprecated name of ``from_array_to_matrix`` (SymPy 1.11)."""
from . import from_array_to_matrix
from .conv_array_to_indexed import _conv_to_from_decorator

convert_array_to_matrix = _conv_to_from_decorator(from_array_to_matrix.convert_array_to_matrix)
_array2matrix = _conv_to_from_decorator(from_array_to_matrix._array2matrix)
_remove_trivial_dims = _conv_to_from_decorator(from_array_to_matrix._remove_trivial_dims)

__all__ = ['convert_array_to_matrix']
