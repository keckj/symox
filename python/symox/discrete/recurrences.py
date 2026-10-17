"""Linear recurrences (``sympy.discrete.recurrences``)."""

from .. import _core
from ..ntheory._intfunc import as_int
from .transforms import iterable


def _ints(seq):
    out = []
    for x in seq:
        if type(x) is int:
            out.append(x)
        elif getattr(x, 'is_Integer', False) is True:
            out.append(int(x))
        else:
            return None
    return out


def linrec(coeffs, init, n):
    """Value at n of the homogeneous linear recurrence with constant
    coefficients ``coeffs`` and initial values ``init``."""
    from .. import S, sympify, Integer
    if not coeffs:
        return S.Zero
    if not iterable(coeffs):
        raise TypeError("Expected a sequence of coefficients for"
                        " the recurrence")
    if not iterable(init):
        raise TypeError("Expected a sequence of values for the "
                        "initialization of the recurrence")
    n = as_int(n)
    if n < 0:
        raise ValueError("Point of evaluation of recurrence must be a "
                         "non-negative integer")
    c = [sympify(arg) for arg in coeffs]
    b = [sympify(arg) for arg in init]
    k = len(c)
    if len(b) > k:
        raise TypeError("Count of initial values should not exceed the "
                        "order of the recurrence")
    b += [S.Zero] * (k - len(b))
    if n < k:
        return b[n]
    ic, ib = _ints(c), _ints(b)
    if ic is not None and ib is not None:
        return Integer(_core.nt_linrec(ic, ib, n))
    terms = [u * v for u, v in zip(linrec_coeffs(c, n), b)]
    return sum(terms[:-1], terms[-1])


def linrec_coeffs(c, n):
    """Coefficients of ``x**n`` reduced by ``x**k = c_0 x**(k-1) + ... +
    c_(k-1)``, as used by ``linrec``."""
    from .. import S, Integer
    ic = _ints(c)
    if ic and n >= 0:
        return [Integer(v) for v in _core.nt_linrec_coeffs(ic, n)]
    k = len(c)

    def _square_and_reduce(u, offset):
        w = [S.Zero] * (2 * len(u) - 1 + offset)
        for i, p in enumerate(u):
            for j, q in enumerate(u):
                w[offset + i + j] += p * q
        for j in range(len(w) - 1, k - 1, -1):
            for i in range(k):
                w[j - i - 1] += w[j] * c[i]
        return w[:k]

    def _final_coeffs(n):
        if n < k:
            return [S.Zero] * n + [S.One] + [S.Zero] * (k - n - 1)
        return _square_and_reduce(_final_coeffs(n // 2), n % 2)

    return _final_coeffs(n)
