"""Hexadecimal digits of pi (``sympy.ntheory.bbp_pi``)."""

from .. import _core
from ._intfunc import as_int


def pi_hex_digits(n, prec=14):
    """``prec`` hexadecimal digits of pi starting at the n-th one."""
    n, prec = as_int(n), as_int(prec)
    if n < 0:
        raise ValueError('n cannot be negative')
    if prec < 0:
        raise ValueError('prec cannot be negative')
    return _core.nt_pi_hex_digits(n, prec)
