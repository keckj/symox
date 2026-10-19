"""Residues modulo integers (``sympy.ntheory.residue_ntheory``).

The algorithms (orders, primitive roots, Tonelli-Shanks, Hensel lifting,
n-th roots, discrete logarithms) run natively and enumerate roots in
SymPy's order.
"""

from .. import _core
from ._intfunc import as_int
import operator


def n_order(a, n):
    """Multiplicative order of a modulo n."""
    a, n = as_int(a), as_int(n)
    return int(_core.nt_n_order(a, n))


def primitive_root(p, smallest=True):
    """A primitive root of p (the smallest by default), or ``None``."""
    p = as_int(p)
    r = _core.nt_primitive_root(p, bool(smallest))
    return None if r is None else int(r)


def is_primitive_root(a, p):
    """Whether a is a primitive root of p."""
    a, p = as_int(a), as_int(p)
    return _core.nt_is_primitive_root(a, p)


def sqrt_mod_iter(a, p, domain=int):
    """Iterate over the solutions of ``x**2 = a mod p``."""
    a, p = as_int(a), abs(as_int(p))
    if p == 0:
        a % p
    for r in _core.nt_sqrt_mod_iter(a, p):
        yield domain(r)


def sqrt_mod(a, p, all_roots=False):
    """A root of ``x**2 = a mod p`` (the smallest), or all sorted roots."""
    if all_roots:
        return sorted(sqrt_mod_iter(a, p))
    p = abs(as_int(p))
    halfp = p // 2
    x = None
    for r in sqrt_mod_iter(a, p):
        if r < halfp:
            return r
        elif r > halfp:
            return p - r
        else:
            x = r
    return x


def is_quad_residue(a, p):
    """Whether a is a quadratic residue modulo p."""
    a, p = as_int(a), as_int(p)
    return _core.nt_is_quad_residue(a, p)


def is_nthpow_residue(a, n, m):
    """Whether ``x**n = a mod m`` has solutions."""
    a = a % m
    a, n, m = as_int(a), as_int(n), as_int(m)
    return _core.nt_is_nthpow_residue(a, n, m)


def nthroot_mod(a, n, p, all_roots=False):
    """Solutions of ``x**n = a mod p`` (the smallest, or all sorted)."""
    a = a % p
    a, n, p = as_int(a), as_int(n), as_int(p)
    if n < 1:
        raise ValueError("n should be positive")
    if p < 1:
        raise ValueError("p should be positive")
    if n == 1:
        return [a] if all_roots else a
    if n == 2:
        return sqrt_mod(a, p, all_roots)
    ret = [int(r) for r in _core.nt_nthroot_mod_all(a, n, p)]
    if all_roots:
        return ret
    if ret:
        return ret[0]
    return None


def quadratic_residues(p):
    """Sorted list of the quadratic residues modulo p."""
    p = as_int(p)
    if p <= 0:
        return sorted({pow(i, 2, p) for i in range(p // 2 + 1)})
    return [int(r) for r in _core.nt_quadratic_residues(p)]


def discrete_log(n, a, b, order=None, prime_order=None):
    """The discrete logarithm of a to the base b modulo n."""
    n, a, b = as_int(n), as_int(a), as_int(b)
    if n < 1:
        raise ValueError("n should be positive")
    if n == 1:
        return 0
    if order is not None:
        order = as_int(order)
    if prime_order is not None:
        prime_order = bool(prime_order)
    return int(_core.nt_discrete_log(n, a, b, order, prime_order))


def quadratic_congruence(a, b, c, n):
    """Sorted solutions of ``a*x**2 + b*x + c = 0 mod n``."""
    a = as_int(a)
    b = as_int(b)
    c = as_int(c)
    n = as_int(n)
    roots, full = _core.nt_quadratic_congruence(a, b, c, n)
    if a % n == 0 and not full:
        from .. import Integer
        return [Integer(r) for r in roots]
    return [int(r) for r in roots]


def _valid_expr(expr):
    """Integer coefficients of a univariate polynomial expression."""
    from .. import Poly, ZZ
    if not expr.is_polynomial():
        raise ValueError("The expression should be a polynomial")
    if len(expr.free_symbols) > 1:
        raise ValueError("The expression should be univariate")
    polynomial = Poly(expr)
    if not polynomial.domain == ZZ:
        raise ValueError("The expression should should have integer "
                         "coefficients")
    return polynomial.all_coeffs()


def gf_csolve(f, n):
    """Sorted roots modulo n of the integer polynomial f (dense, leading
    coefficient first)."""
    return [int(r) for r in _core.nt_csolve([int(c) for c in f], int(n))]


def polynomial_congruence(expr, m):
    """Sorted solutions of the polynomial congruence ``expr = 0 mod m``."""
    coefficients = _valid_expr(expr)
    from .. import Integer
    coefficients = [Integer(int(num) % int(m)) for num in coefficients]
    rank = len(coefficients)
    if rank == 3:
        return quadratic_congruence(*coefficients, m)
    if rank == 2:
        return quadratic_congruence(0, *coefficients, m)
    if coefficients[0] == 1 and 1 + coefficients[-1] == sum(coefficients):
        return nthroot_mod(-int(coefficients[-1]), rank - 1, int(m), True)
    return gf_csolve(coefficients, m)


def binomial_mod(n, m, k):
    """``binomial(n, m) % k`` (Lucas, Granville and the CRT)."""
    if k < 1:
        raise ValueError('k is required to be positive')
    if n < 0 or m < 0 or m > n:return 0
    k=as_int(k)
    if k==1:return 0
    # Do not silently truncate nonintegral data into a different problem.
    return int(_core.nt_binomial_mod(operator.index(n),operator.index(m),k))


def _head(name):
    def f(*args):
        import symox
        return getattr(symox, name)(*args)
    f.__name__ = name
    f.__doc__ = "Deprecated in SymPy: the ``%s`` function head." % name
    return f


legendre_symbol = _head('legendre_symbol')
jacobi_symbol = _head('jacobi_symbol')
mobius = _head('mobius')
