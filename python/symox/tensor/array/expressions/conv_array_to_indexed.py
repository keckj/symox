"""Deprecated name of ``from_array_to_indexed`` (SymPy 1.11): the functions
warn, then call the ``from_`` module."""
from . import from_array_to_indexed
from ....utilities.decorator import deprecated


_conv_to_from_decorator = deprecated(
    "module has been renamed by replacing 'conv_' with 'from_' in its name",
    deprecated_since_version="1.11",
    active_deprecations_target="deprecated-conv-array-expr-module-names",
)


convert_array_to_indexed = _conv_to_from_decorator(from_array_to_indexed.convert_array_to_indexed)

__all__ = ['convert_array_to_indexed']
