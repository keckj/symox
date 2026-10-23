"""mpmath's ``ellippi`` with a controlled initial integral.

SymPy 1.14 evaluates ``elliptic_pi`` with ``mpmath.ellippi`` at the target
precision. For a complete or incomplete integral whose Carlson ``R_J``
arguments leave the region where Carlson's duplication is valid (``m > 1``,
``n > 1`` and most complex arguments), mpmath's ``R_J`` first integrates
numerically along the segment ``[0, N + i*margin]`` with ``quadsubdiv`` and
then applies Carlson's algorithm to the shifted arguments. That quadrature
passes close to the integrand's singularities: it can take minutes at some
precisions and it returns wrong trailing digits at others.

:func:`ellippi` follows mpmath's ``ellippi`` and ``elliprj`` step by step,
with the same integration path, the same principal branches, the same
special values and the same Carlson tail. Only the segment integral differs:
it is split until every piece is far from the singular points ``-x, -y, -z,
-p`` (with the endpoint singularity at ``t = 0`` removed by ``t = b*u**2``),
and each piece is integrated by Gauss-Legendre with an error check. The value
is therefore the same analytic value as mpmath's, to the working precision.
"""

import cmath

import mpmath.functions.elliptic as _elliptic
from mpmath import mp
from mpmath.calculus.quadrature import GaussLegendre

__all__ = ["ellippi", "elliprj"]

# Maximum bisection depth of a segment; beyond it the integral is left to
# mpmath's own quadrature.
_MAX_DEPTH = 64


def _distance(q, a, b):
    """Distance from the point ``q`` to the segment ``[a, b]``."""
    d = b - a
    s = ((q - a) * mp.conj(d)).real / abs(d) ** 2
    s = min(max(s, 0), 1)
    return abs(q - (a + s * d))


class _Unresolved(Exception):
    pass


# Gauss-Legendre rules on [-1, 1] with 3*2**(k-1) nodes, kept at the highest
# precision computed so far (at least 384 bits, in steps of 128): a rule
# computed at more bits serves any lower precision, so the successive
# precisions of one evaluation share it. Above 48 nodes, bisecting the piece
# is cheaper than computing the next rule.
_RULES = {}
_MAX_LEVEL = 5


def _rule(level, prec):
    cached = _RULES.get(level)
    if cached is None or cached[0] < prec:
        bits = max(384, -(-prec // 128) * 128)
        orig = mp.prec
        try:
            mp.prec = bits + 20
            nodes = GaussLegendre(mp).calc_nodes(level, bits)
        finally:
            mp.prec = orig
        cached = _RULES[level] = (bits, nodes)
    return cached[1]


def _gauss_legendre(f, a, b, tol):
    """Integral of ``f`` on ``[a, b]``, or ``None`` if not converged.

    On a piece four times its length away from the singular points, the
    error of the rule with ``n`` nodes decreases like ``rho**(-2*n)``:
    doubling the nodes squares the relative error. Once the differences of
    successive rules decrease that way, the difference ``d`` of the last two
    rules bounds the error of the previous one, and the last one is accurate
    to about ``d**2``; it is accepted when ``d**2`` is 1024 times below the
    tolerance.
    """
    half = (b - a) / 2
    mid = (a + b) / 2
    previous = None
    last = None
    for level in range(1, _MAX_LEVEL + 1):
        value = half * mp.fsum(w * f(mid + half * x) for x, w in _rule(level, mp.prec))
        if previous is not None:
            size = abs(value)
            difference = abs(value - previous)
            if difference <= tol * size:
                return value
            relative = difference / size if size else mp.inf
            if last is not None and relative <= last**1.5 and relative**2 * 1024 <= tol:
                return value
            last = relative
        previous = value
    return None


def _integral(x, y, z, p, c):
    """``1.5 * integral(F, 0, c)`` along the segment ``[0, c]``.

    ``F(t) = 1/(sqrt(t+x)*sqrt(t+y)*sqrt(t+z)*(t+p))`` with principal square
    roots, as in mpmath's ``RJ_calc``. The segment is bisected until every
    piece is at least four times its length away from the singular points, then
    each piece is integrated by Gauss-Legendre and bisected again while the
    error estimate exceeds the working precision. A zero among ``x, y, z``
    (at most one) puts a ``1/sqrt(t)`` singularity at ``t = 0``, removed on
    the first piece by ``t = b*u**2``.
    """
    sqrt = mp.sqrt
    ts = [x, y, z]
    zero = any(not w for w in ts)
    # The zero point first when one of x, y, z vanishes.
    singular = ([mp.zero] if zero else []) + [-w for w in ts + [p] if w]

    def F(t):
        a, b, c = t + x, t + y, t + z
        # One square root instead of three: sqrt(a)*sqrt(b)*sqrt(c) is
        # +-sqrt(a*b*c), the sign read off a double-precision product (the
        # two candidates are opposite, far apart at any precision).
        root = sqrt(a * b * c)
        try:
            low = cmath.sqrt(complex(a)) * cmath.sqrt(complex(b)) * cmath.sqrt(complex(c))
        except (OverflowError, ValueError):
            low = None
        if low is None or not cmath.isfinite(low) or not low:
            return 1 / (sqrt(a) * sqrt(b) * sqrt(c) * (t + p))
        if (low * complex(root).conjugate()).real < 0:
            root = -root
        return 1 / (root * (t + p))

    def G(b, u):
        # sqrt(b*u**2) = u*sqrt(b) for u > 0.
        t = b * u * u
        rest = [sqrt(t + w) for w in ts if w]
        return 2 * sqrt(b) / (rest[0] * rest[1] * (t + p))

    tol = mp.eps * 16

    def piece(a, b, depth):
        # The zero singular point is removed on the first piece only, where
        # u = sqrt(t/b) brings the other ones closer: four times as far in u
        # needs sixteen times as far in t.
        first = not a and zero
        points = singular[1:] if first else singular
        if min(_distance(q, a, b) for q in points) >= (16 if first else 4) * abs(b - a):
            if first:
                value = _gauss_legendre(lambda u: G(b, u), mp.zero, mp.one, tol)
            else:
                value = _gauss_legendre(F, a, b, tol)
            if value is not None:
                return value
        if depth >= _MAX_DEPTH:
            raise _Unresolved
        m = (a + b) / 2
        return piece(a, m, depth + 1) + piece(m, b, depth + 1)

    return 1.5 * piece(mp.zero, c, 0)


def _integrated(x, y, z, p, r):
    """mpmath's integration branch of ``RJ_calc``, or ``None``.

    mpmath integrates along ``[0, c]``, ``c = N + margin``, and applies
    Carlson's algorithm to the arguments shifted by ``c``. When ``margin`` is
    ``i`` (or ``-i``), every singular point ``-w`` lies below (above) that
    path or to the left of the imaginary axis, and no branch cut of
    ``sqrt(t + w)`` (the ray ``-w - [0, oo)``) meets the region between the
    path and the half-line ``i*sigma*[0, oo)``. The integral is then taken
    on that half-line instead: along ``[0, i*sigma*H]`` numerically, and
    beyond it, with ``t = i*sigma*(H + u)``, as
    ``exp(-3*i*sigma*pi/4) * R_J(H - i*sigma*w)``, where every argument has
    a positive real part, so that Carlson's algorithm is valid.
    """
    ts = [x, y, z, p]
    N = mp.ceil(-min(t.real for t in ts)) + 1
    if all((t.imag >= 0 or t.real > 0) for t in ts):
        sigma = 1
    elif all((t.imag < 0 or t.real > 0) for t in ts):
        sigma = -1
    else:
        # mpmath's path between singular points of both half-planes.
        margin = 1
        for t in ts:
            if t.imag >= 0 or t.real > 0:
                continue
            margin = min(margin, abs(t.imag) * 0.5)
        c = N + margin * mp.j
        initial = _integral(x, y, z, p, c)
        return initial + _elliptic.RJ_calc(mp, x + c, y + c, z + c, p + c, r, 0)
    rotation = sigma * mp.j
    low = max(-sigma * t.imag for t in ts)
    sizes = [abs(t) for t in ts if t]
    H = max(low, 0) + min(1, min(sizes) / 4)
    initial = _integral(x, y, z, p, rotation * H)
    shifted = [H - rotation * t for t in ts]
    factor = mp.expjpi(-sigma * mp.mpf(3) / 4)
    return initial + factor * _elliptic.RJ_calc(mp, *shifted, r, 0)


def _rj_calc(x, y, z, p, r):
    """mpmath's ``RJ_calc(..., integration=1)``, initial integral replaced."""
    ctx = mp
    if not (ctx.isnormal(x) and ctx.isnormal(y) and ctx.isnormal(z) and ctx.isnormal(p)):
        if ctx.isnan(x) or ctx.isnan(y) or ctx.isnan(z) or ctx.isnan(p):
            return x * y * z
        if ctx.isinf(x) or ctx.isinf(y) or ctx.isinf(z) or ctx.isinf(p):
            return ctx.zero
    if not p:
        return ctx.inf
    if (not x) + (not y) + (not z) > 1:
        return ctx.inf
    ok = (x.real >= 0 and y.real >= 0 and z.real >= 0 and p.real > 0)
    if not ok:
        if x == p or y == p or z == p:
            ok = True
    if not ok:
        if p.imag != 0 or p.real >= 0:
            if (x.imag == 0 and x.real >= 0 and ctx.conj(y) == z):
                ok = True
            if (y.imag == 0 and y.real >= 0 and ctx.conj(x) == z):
                ok = True
            if (z.imag == 0 and z.real >= 0 and ctx.conj(x) == y):
                ok = True
    if ok:
        return _elliptic.RJ_calc(ctx, x, y, z, p, r, 1)
    try:
        return _integrated(x, y, z, p, r)
    except _Unresolved:
        # A singular point on the path: mpmath's own quadrature.
        return _elliptic.RJ_calc(ctx, x, y, z, p, r, 1)


def elliprj(x, y, z, p):
    """mpmath's ``elliprj`` (``integration=1``)."""
    ctx = mp
    x = ctx.convert(x)
    y = ctx.convert(y)
    z = ctx.convert(z)
    p = ctx.convert(p)
    prec = ctx.prec
    try:
        ctx.prec += 20
        tol = ctx.eps * 2**10
        v = _rj_calc(x, y, z, p, tol)
    finally:
        ctx.prec = prec
    return +v


def ellippi(*args):
    """mpmath's ``ellippi`` through :func:`elliprj`."""
    ctx = mp
    args = [ctx.convert(arg) for arg in args]
    prec = ctx.prec
    try:
        # mpmath's defun_wrapped: ten guard bits, rounded on return.
        ctx.prec += 10
        value = _ellippi(ctx, *args)
    finally:
        ctx.prec = prec
    return +value


def _ellippi(ctx, *args):
    if len(args) == 2:
        n, m = args
        complete = True
        z = phi = ctx.pi / 2
    else:
        n, phi, m = args
        complete = False
        z = phi
    if not (ctx.isnormal(n) and ctx.isnormal(z) and ctx.isnormal(m)):
        if ctx.isnan(n) or ctx.isnan(z) or ctx.isnan(m):
            raise ValueError
        if complete:
            if m == 0:
                if n == 1:
                    return ctx.inf
                return ctx.pi / (2 * ctx.sqrt(1 - n))
            if n == 0:
                return ctx.ellipk(m)
            if ctx.isinf(n) or ctx.isinf(m):
                return ctx.zero
        else:
            if z == 0:
                return z
            if ctx.isinf(n):
                return ctx.zero
            if ctx.isinf(m):
                return ctx.zero
        if ctx.isinf(n) or ctx.isinf(z) or ctx.isinf(m):
            raise ValueError
    if complete:
        if m == 1:
            if n == 1:
                return ctx.inf
            return -ctx.inf / ctx.sign(n - 1)
        away = False
    else:
        x = z.real
        ctx.prec += max(0, ctx.mag(x))
        pi = +ctx.pi
        away = abs(x) > pi / 2
    if away:
        d = ctx.nint(x / pi)
        z = z - pi * d
        P = 2 * d * ellippi(n, m)
        if ctx.isinf(P):
            return ctx.inf
    else:
        P = 0

    def terms():
        if complete:
            c, s = ctx.zero, ctx.one
        else:
            c, s = ctx.cos_sin(z)
        x = c**2
        y = 1 - m * s**2
        RF = ctx.elliprf(x, y, 1)
        RJ = elliprj(x, y, 1, 1 - n * s**2)
        return s * RF, n * s**3 * RJ / 3

    return ctx.sum_accurately(terms) + P
