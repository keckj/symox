"""Shape queries for matrices, arrays and user registered expression types."""
from functools import singledispatch


class NoShapeError(ValueError):
    """An expression has no matrix or array shape."""


@singledispatch
def shape(expr):
    """Return an object's shape; additional types can use ``shape.register``."""
    if hasattr(expr, 'shape'):
        return expr.shape
    raise NoShapeError('%s does not have shape, or its type is not registered to shape().' % expr)


__all__ = ['shape']
