"""``sympy.matrices.utilities``: the ``dotprodsimp`` switch and the zero
tests and simplification wrappers shared by the matrix routines."""
from contextlib import contextmanager
from threading import local


class DotProdSimpState(local):
    def __init__(self):
        self.state = None


_dotprodsimp_state = DotProdSimpState()


@contextmanager
def dotprodsimp(x):
    old = _dotprodsimp_state.state

    try:
        _dotprodsimp_state.state = x
        yield
    finally:
        _dotprodsimp_state.state = old


def _dotprodsimp(expr, withsimp=False):
    from symox.simplify.simplify import dotprodsimp as dps
    return dps(expr, withsimp=withsimp)


def _get_intermediate_simp(deffunc=lambda x: x, offfunc=lambda x: x,
        onfunc=_dotprodsimp, dotprodsimp=None):
    if dotprodsimp is False or _dotprodsimp_state.state is False:
        return offfunc
    if dotprodsimp is True or _dotprodsimp_state.state is True:
        return onfunc

    return deffunc


def _get_intermediate_simp_bool(default=False, dotprodsimp=None):
    return _get_intermediate_simp(default, False, True, dotprodsimp)


def _iszero(x):
    """Returns True if x is zero."""
    return getattr(x, 'is_zero', None)


def _is_zero_after_expand_mul(x):
    """Tests by expand_mul only, suitable for polynomials and rational
    functions."""
    from symox import expand_mul
    return expand_mul(x) == 0


def _simplify(expr):
    from symox import simplify
    return simplify(expr)
