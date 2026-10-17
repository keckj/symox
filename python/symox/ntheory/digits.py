"""Digit expansions (``sympy.ntheory.digits``)."""

from collections import defaultdict

from ._intfunc import as_int


def digits(n, b=10, digits=None):
    """``[b, d0, d1, ...]`` (``-b`` first when n is negative)."""
    b = as_int(b)
    n = as_int(n)
    if b < 2:
        raise ValueError("b must be greater than 1")
    x, y = abs(n), []
    while x >= b:
        x, r = divmod(x, b)
        y.append(r)
    y.append(x)
    y.append(-b if n < 0 else b)
    y.reverse()
    ndig = len(y) - 1
    if digits is not None:
        if ndig > digits:
            raise ValueError(
                "For %s, at least %s digits are needed." % (n, ndig))
        elif ndig < digits:
            y[1:1] = [0] * (digits - ndig)
    return y


def count_digits(n, b=10):
    """A ``defaultdict(int)`` mapping each digit of n to its count."""
    counts = {}
    for d in digits(n, b):
        counts[d] = counts.get(d, 0) + 1
    rv = defaultdict(int, counts.items())
    rv.pop(b) if b in rv else rv.pop(-b)
    return rv


def is_palindromic(n, b=10):
    """Whether the base-b digits of n read the same in both directions."""
    s = digits(n, b)[1:]
    return s == s[::-1]
