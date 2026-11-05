"""Deprecated name of ``from_matrix_to_array`` (SymPy 1.11)."""
from . import from_matrix_to_array
from .conv_array_to_indexed import _conv_to_from_decorator

convert_matrix_to_array = _conv_to_from_decorator(from_matrix_to_array.convert_matrix_to_array)

__all__ = ['convert_matrix_to_array']
