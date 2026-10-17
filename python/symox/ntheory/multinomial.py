"""Binomial and multinomial coefficients (``sympy.ntheory.multinomial``)."""

from .. import _core
from ._intfunc import as_int


def binomial_coefficients(n):
    """Dict of ``(k1, k2): C(n, k1)`` with ``k1 + k2 == n``."""
    n = as_int(n)
    d = {(0, n): 1, (n, 0): 1}
    if n > 1:
        row = _core.nt_binomial_row(n)
        for k in range(1, n // 2 + 1):
            d[k, n - k] = d[n - k, k] = int(row[k])
    return d


def binomial_coefficients_list(n):
    """Row ``n`` of Pascal's triangle as a list."""
    n = as_int(n)
    if n < 0:
        return [1] * (n + 1)
    return [int(c) for c in _core.nt_binomial_row(n)]


def multinomial_coefficients(m, n):
    """Dict of multinomial coefficients ``(k1, ..., km): n!/(k1!...km!)``."""
    m = as_int(m)
    n = as_int(n)
    if not m:
        if n:
            return {}
        return {(): 1}
    if m == 2:
        return binomial_coefficients(n)
    if m >= 2 * n and n > 1:
        return dict(multinomial_coefficients_iterator(m, n))
    t = [n] + [0] * (m - 1)
    r = {tuple(t): 1}
    j = 0 if n else m
    while j < m - 1:
        tj = t[j]
        if j:
            t[j] = 0
            t[0] = tj
        if tj > 1:
            t[j + 1] += 1
            j = 0
            start = 1
            v = 0
        else:
            j += 1
            start = j + 1
            v = r[tuple(t)]
            t[j] += 1
        for k in range(start, m):
            if t[k]:
                t[k] -= 1
                v += r[tuple(t)]
                t[k] += 1
        t[0] -= 1
        r[tuple(t)] = v * tj // (n - t[0])
    return r


def multinomial_coefficients_iterator(m, n, _tuple=tuple):
    """Iterator over ``multinomial_coefficients(m, n)`` items."""
    m = as_int(m)
    n = as_int(n)
    if m < 2 * n or n == 1:
        mc = multinomial_coefficients(m, n)
        yield from mc.items()
        return
    mc = multinomial_coefficients(n, n)
    mc = {_tuple(filter(None, k)): v for k, v in mc.items()}
    t = [n] + [0] * (m - 1)
    t1 = _tuple(t)
    yield (t1, mc[_tuple(filter(None, t1))])
    j = 0 if n else m
    while j < m - 1:
        tj = t[j]
        if j:
            t[j] = 0
            t[0] = tj
        if tj > 1:
            t[j + 1] += 1
            j = 0
        else:
            j += 1
            t[j] += 1
        t[0] -= 1
        t1 = _tuple(t)
        yield (t1, mc[_tuple(filter(None, t1))])
