"""Deprecated name of ``from_indexed_to_array`` (SymPy 1.11)."""
from . import from_indexed_to_array
from .conv_array_to_indexed import _conv_to_from_decorator

convert_indexed_to_array = _conv_to_from_decorator(from_indexed_to_array.convert_indexed_to_array)

__all__ = ['convert_indexed_to_array']
