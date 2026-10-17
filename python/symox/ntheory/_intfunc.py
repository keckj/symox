"""Integer helpers matching ``sympy.core.intfunc`` and ``as_int``.

The arithmetic runs on Python integers or in the native core; this module
only reproduces SymPy's argument conventions and error messages.
"""

import math
import operator

from .. import _core


def as_int(n, strict=True):
    """``sympy.utilities.misc.as_int``."""
    if strict:
        try:
            if isinstance(n, bool):
                raise TypeError
            return operator.index(n)
        except TypeError:
            raise ValueError('%s is not an integer' % (n,))
    try:
        result = int(n)
    except TypeError:
        raise ValueError('%s is not an integer' % (n,))
    if n - result:
        raise ValueError('%s is not an integer' % (n,))
    return result


def _int_ceiling(n):
    """``int(ceiling(n))`` for numbers (``_as_int_ceiling``)."""
    try:
        return operator.index(n)
    except TypeError:
        pass
    from .. import ceiling, sympify
    return int(ceiling(sympify(n)))


def _int_floor(n):
    try:
        return operator.index(n)
    except TypeError:
        return int(n)


def _bit_scan1(n):
    n = abs(n)
    return (n & -n).bit_length() - 1


def trailing(n):
    """Number of trailing zero bits of ``n``."""
    if not n:
        return 0
    return _bit_scan1(int(n))


def integer_log(n, b):
    """``(e, bool)`` with ``e`` the largest integer such that ``|b**e| <= |n|``."""
    n = as_int(n)
    b = as_int(b)
    if b < 0:
        e, t = integer_log(abs(n), -b)
        t = t and e % 2 == (n < 0)
        return e, t
    if b <= 1:
        raise ValueError('base must be 2 or more')
    if n < 0:
        if b != 2:
            return 0, False
        e, t = integer_log(-n, b)
        return e, False
    if n == 0:
        raise ValueError('n cannot be 0')
    return _core.nt_integer_log(n, b)


def num_digits(n, base=10):
    """Number of digits of ``n`` in ``base``."""
    if base < 0:
        raise ValueError('base must be int greater than 1')
    if not n:
        return 1
    e, t = integer_log(abs(n), base)
    return 1 + e


def igcd(*args):
    """Integer gcd (as a Python ``int``)."""
    if len(args) < 2:
        raise TypeError('igcd() takes at least 2 arguments (%s given)' % len(args))
    return int(math.gcd(*map(as_int, args)))


def ilcm(*args):
    """Integer lcm (as a Python ``int``)."""
    if len(args) < 2:
        raise TypeError('ilcm() takes at least 2 arguments (%s given)' % len(args))
    return int(math.lcm(*map(as_int, args)))


igcd2 = math.gcd


def igcdex(a, b):
    """``(x, y, g)`` with ``x*a + y*b == g``."""
    return tuple(_core.nt_igcdex(int(a), int(b)))


def mod_inverse(a, m):
    """Inverse of ``a`` modulo ``m`` (``1/a`` for numbers when ``m > 1``)."""
    c = None
    try:
        a, m = as_int(a), as_int(m)
        if m != 1 and m != -1:
            x, _, g = igcdex(a, m)
            if g == 1:
                c = x % m
    except ValueError:
        from .. import S, sympify
        a, m = sympify(a), sympify(m)
        if not (a.is_number and m.is_number):
            raise TypeError(
                "Expected numbers for arguments; symbolic `mod_inverse` is not\n"
                "implemented but symbolic expressions can be handled with the\n"
                "similar function, sympy.polys.polytools.invert")
        big = m > 1
        if big not in (S.true, S.false, True, False):
            raise ValueError('m > 1 did not evaluate; try to simplify %s' % m)
        elif big:
            c = 1 / a
    if c is None:
        raise ValueError('inverse of %s (mod %s) does not exist' % (a, m))
    return c


def integer_nthroot(y, n):
    """``(floor(y**(1/n)), exact)``."""
    y, n = as_int(y), as_int(n)
    if y < 0:
        raise ValueError("y must be nonnegative")
    if n < 1:
        raise ValueError("n must be positive")
    x, b = _core.nt_iroot(y, n)
    return int(x), b


def isqrt(n):
    if n < 0:
        raise ValueError('n must be nonnegative')
    return math.isqrt(int(n))
