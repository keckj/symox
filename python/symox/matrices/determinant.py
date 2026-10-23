"""Pivot searches of ``sympy.matrices.determinant``."""
from .._matrix_definite import _find_reasonable_pivot
from .utilities import _iszero

__all__ = ['_find_reasonable_pivot', '_find_reasonable_pivot_naive']


def _find_reasonable_pivot_naive(col, iszerofunc=_iszero, simpfunc=None):
    """Find the lowest index of an item in ``col`` that is suitable for a
    pivot.  Returns ``(index, value, assumed_nonzero, newly_determined)``
    exactly as SymPy's ``_find_reasonable_pivot_naive``."""
    indeterminates = []
    for i, col_val in enumerate(col):
        col_val_is_zero = iszerofunc(col_val)
        if col_val_is_zero == False:  # noqa: E712
            return i, col_val, False, []
        elif col_val_is_zero is None:
            indeterminates.append((i, col_val))

    if len(indeterminates) == 0:
        return None, None, False, []

    if simpfunc is None:
        return indeterminates[0][0], indeterminates[0][1], True, []

    newly_determined = []
    for i, col_val in indeterminates:
        tmp_col_val = simpfunc(col_val)
        if id(col_val) != id(tmp_col_val):
            newly_determined.append((i, tmp_col_val))
            if iszerofunc(tmp_col_val) == False:  # noqa: E712
                return i, tmp_col_val, False, newly_determined

    return indeterminates[0][0], indeterminates[0][1], True, newly_determined
