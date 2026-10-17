"""Orthogonal polynomials (``sympy.polys.orthopolys``).

The three-term recurrences run natively (``symox::poly::orthopolys``) over
SymPy's ground domain: ``ZZ``/``QQ`` for the parameter-free families, ``QQ``,
``RR`` or the fraction field ``ZZ(a, b)`` of the parameters otherwise. Other
parameters (algebraic or transcendental constants, floats mixed with symbols)
use canceled exact expressions and an inferred domain.
"""
from .._core import Poly, _orthopoly_coefficients
from .specialpolys import _named_result

__all__ = [
    'chebyshevt_poly', 'chebyshevu_poly', 'gegenbauer_poly', 'hermite_poly',
    'hermite_prob_poly', 'jacobi_poly', 'laguerre_poly', 'legendre_poly',
    'spherical_bessel_fn',
]


def _named(family, n, name, x, params, polys):
    """SymPy's ``named_poly``: validate ``n`` then build the result."""
    if n < 0:
        raise ValueError("Cannot generate %s of index %s" % (name, n))
    coeffs, domain = _orthopoly_coefficients(family, int(n), list(params))
    if not coeffs:
        from .. import S
        coeffs = [S.Zero]
    return _named_result(coeffs, x, domain, polys)


def jacobi_poly(n, a, b, x=None, polys=False):
    """Jacobi polynomial :math:`P_n^{(a,b)}(x)`."""
    return _named('jacobi', n, "Jacobi polynomial", x, (a, b), polys)


def gegenbauer_poly(n, a, x=None, polys=False):
    """Gegenbauer polynomial :math:`C_n^{(a)}(x)`."""
    return _named('gegenbauer', n, "Gegenbauer polynomial", x, (a,), polys)


def chebyshevt_poly(n, x=None, polys=False):
    """Chebyshev polynomial of the first kind :math:`T_n(x)`."""
    return _named('chebyshevt', n, "Chebyshev polynomial of the first kind", x, (), polys)


def chebyshevu_poly(n, x=None, polys=False):
    """Chebyshev polynomial of the second kind :math:`U_n(x)`."""
    return _named('chebyshevu', n, "Chebyshev polynomial of the second kind", x, (), polys)


def hermite_poly(n, x=None, polys=False):
    """Physicist's Hermite polynomial :math:`H_n(x)`."""
    return _named('hermite', n, "Hermite polynomial", x, (), polys)


def hermite_prob_poly(n, x=None, polys=False):
    """Probabilist's Hermite polynomial :math:`He_n(x)`."""
    return _named('hermite_prob', n, "probabilist's Hermite polynomial", x, (), polys)


def legendre_poly(n, x=None, polys=False):
    """Legendre polynomial :math:`P_n(x)`."""
    return _named('legendre', n, "Legendre polynomial", x, (), polys)


def laguerre_poly(n, x=None, alpha=0, polys=False):
    """Generalized Laguerre polynomial :math:`L_n^{(alpha)}(x)`."""
    return _named('laguerre', n, "Laguerre polynomial", x, (alpha,), polys)


def spherical_bessel_fn(n, x=None, polys=False):
    """Coefficients ``fn(n, z)`` of the spherical Bessel functions, a
    polynomial in ``1/z`` with ``fn(n-1) + fn(n+1) = (2*n+1)/z * fn(n)``."""
    from .. import Dummy, Integer
    if x is None:
        x = Dummy('x')
    family = 'spherical_bessel_fn_minus' if n < 0 else 'spherical_bessel_fn'
    coeffs, domain = _orthopoly_coefficients(family, int(abs(n)), [])
    if not coeffs:
        coeffs = [Integer(0)]
    return _named_result(coeffs, Integer(1)/x, domain, polys)
