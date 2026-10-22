"""Shape queries for matrices, arrays and user registered expression types."""
from functools import singledispatch


class NoShapeError(Exception):
    """An expression has no matrix or array shape."""


@singledispatch
def shape(expr):
    """Return an object's shape; additional types can use ``shape.register``."""
    if hasattr(expr, 'shape'):
        return expr.shape
    raise NoShapeError('%s does not have shape, or its type is not registered to shape().' % expr)


def _register_integral():
    from ..integrals import Integral

    @shape.register(Integral)
    def _(expr):
        return shape(expr.function)


_register_integral()


__all__ = ['shape']
