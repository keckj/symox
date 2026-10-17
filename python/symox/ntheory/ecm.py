"""Lenstra's elliptic curve factorization (``sympy.ntheory.ecm``)."""

from .. import _core
from ._intfunc import as_int


class Point:
    """Montgomery curve point ``(X : Z)`` with ``a_24 = (a + 2)/4`` mod ``mod``."""

    def __init__(self, x_cord, z_cord, a_24, mod):
        self.x_cord = x_cord
        self.z_cord = z_cord
        self.a_24 = a_24
        self.mod = mod

    def _t(self):
        return (self.x_cord, self.z_cord, self.a_24, self.mod)

    def __eq__(self, other):
        if self.a_24 != other.a_24 or self.mod != other.mod:
            return False
        return self.x_cord * other.z_cord % self.mod == \
            other.x_cord * self.z_cord % self.mod

    __hash__ = None

    def add(self, Q, diff):
        """Differential addition ``self + Q`` knowing ``diff = self - Q``."""
        return Point(*_core.nt_point_add(self._t(), Q._t(), diff._t()))

    def double(self):
        """``2*self``."""
        return Point(*_core.nt_point_double(self._t()))

    def mont_ladder(self, k):
        """``k*self`` by the Montgomery ladder."""
        return Point(*_core.nt_point_ladder(self._t(), as_int(k)))


def ecm(n, B1=10000, B2=100000, max_curve=200, seed=1234):
    """Set of prime factors of n found by trial division and ECM."""
    n = as_int(n)
    return set(_core.nt_ecm(n, as_int(B1), as_int(B2), as_int(max_curve),
                            seed if seed is not None else 0))
