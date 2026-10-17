"""Egyptian fractions (``sympy.ntheory.egyptian_fraction``)."""

from .. import _core


def egyptian_fraction(r, algorithm='Greedy'):
    """Denominators of an Egyptian fraction expansion of the rational r."""
    from .. import Integer, Rational, Tuple
    if not getattr(r, 'is_Rational', False) or isinstance(r, (int, float)):
        if isinstance(r, (Tuple, tuple)) and len(r) == 2:
            r = Rational(*r)
        else:
            raise ValueError('Value must be a Rational or tuple of ints')
    if r <= 0:
        raise ValueError('Value must be positive')
    x, y = r.as_numer_denom()
    out = _core.nt_egyptian_fraction(int(x), int(y), algorithm)
    return [Integer(i) for i in out]


def egypt_greedy(x, y):
    """Greedy (Fibonacci-Sylvester) expansion of x/y, ``gcd(x, y) == 1``."""
    return [int(i) for i in _core.nt_egypt("Greedy", int(x), int(y))]


def egypt_graham_jewett(x, y):
    """Graham-Jewett expansion of x/y."""
    return [int(i) for i in _core.nt_egypt("Graham Jewett", int(x), int(y))]


def egypt_takenouchi(x, y):
    """Takenouchi expansion of x/y."""
    return [int(i) for i in _core.nt_egypt("Takenouchi", int(x), int(y))]


def egypt_golomb(x, y):
    """Golomb expansion of x/y."""
    return [int(i) for i in _core.nt_egypt("Golomb", int(x), int(y))]


def egypt_harmonic(r):
    """``(prefix, remainder)``: the longest harmonic prefix ``1, 2, ..., d``
    whose sum stays at most r, and ``r`` minus that sum."""
    from .. import Integer, Rational
    prefix, num, den = _core.nt_egypt_harmonic(int(r.p), int(r.q))
    return ([Integer(i) for i in prefix], Rational(num, den))
