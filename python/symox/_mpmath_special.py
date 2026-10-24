"""Multiprecision evaluation of special functions through mpmath.

The native evaluator calls :func:`evaluate` for a special function head that
has no native multiprecision kernel. It mirrors SymPy 1.14's
``Function._eval_evalf`` and the per-class ``_eval_evalf`` overrides: the
same mpmath function, the same argument precision and the same working
precision. Arguments arrive as exact binary values ``(re_man, re_exp,
im_man, im_exp)`` (or nested lists of them for tuple parameters), and the
result is returned in the same form. ``None`` means the value is not
available, and the expression stays unevaluated; nothing is computed in
double precision.
"""

import mpmath
from mpmath import mp, workprec
from mpmath.libmp import from_man_exp

from ._ellippi import ellippi as _ellippi

__all__ = ["evaluate"]


def _real(man, exp):
    return mp.make_mpf(from_man_exp(man, exp))


def _number(value):
    re_man, re_exp, im_man, im_exp = value
    re = _real(re_man, re_exp)
    if not im_man:
        return re
    return mp.make_mpc((re._mpf_, _real(im_man, im_exp)._mpf_))


def _convert(value):
    if isinstance(value, list):
        return [_convert(item) for item in value]
    return _number(value)


def _parts(raw):
    sign, man, exp, bc = raw
    if not man:
        # Zero is (0, 0, 0, 0); infinities and nan have a zero mantissa
        # with a nonzero bit count.
        return (0, 0) if bc == 0 and exp == 0 else None
    return (-man if sign else man, exp)


def _result(value):
    if isinstance(value, mpmath.mpf):
        parts = _parts(value._mpf_)
        return None if parts is None else parts + (0, 0)
    if isinstance(value, mpmath.mpc):
        re, im = (_parts(part) for part in value._mpc_)
        return None if re is None or im is None else re + im
    return None


def _isint(value):
    return mp.isint(value)


def _polygamma(s, z):
    # SymPy's polygamma._eval_evalf, at prec + 12. At a pole it returns zoo,
    # whose evalf is nan (the evaluation fails and stays unevaluated; only
    # an exact Float argument gives nan).
    if _isint(z) and z <= 0:
        return mp.nan
    if _isint(s) and s >= 0:
        return mp.polygamma(s, z)
    zt = mp.zeta(s + 1, z)
    dzt = mp.zeta(s + 1, z, 1)
    return (dzt + (mp.euler + mp.digamma(-s)) * zt) * mp.rgamma(-s)


def _bernoulli(n, x=None):
    x = mp.mpf(1) if x is None else x
    if n == 0:
        return mp.mpf(1)
    if n == 1:
        return x - mp.mpf(0.5)
    if _isint(n) and n >= 0:
        return mp.bernoulli(n) if x == 1 else mp.bernpoly(n, x)
    return -n * mp.zeta(1 - n, x)


def _harmonic(n, m=None):
    m = mp.mpf(1) if m is None else m
    if _isint(n) and n < 0:
        return None
    if m == 1:
        return mp.harmonic(n)
    return mp.zeta(m) - mp.zeta(m, n + 1)


def _euler(m, x=None):
    if _isint(m) and m >= 0:
        return mp.eulernum(m) if x is None else mp.eulerpoly(m, x)
    if m == -1:
        result = mp.pi if x is None else mp.digamma((x + 1) / 2) - mp.digamma(x / 2)
    else:
        y = mp.mpf(0.5) if x is None else x
        result = 2 * (mp.zeta(-m, y) - 2 ** (m + 1) * mp.zeta(-m, (y + 1) / 2))
    if x is None:
        result *= 2 ** m
    return result


def _andre(s):
    sp, cp = mp.sinpi(s / 2), mp.cospi(s / 2)
    return 2 * mp.dirichlet(-s, (-sp, cp, sp, -cp))


def _spherical(kind, n, z):
    # jn/yn rewrite to besselj/bessely for an integer order.
    if not _isint(n):
        return None
    order = n + mp.mpf(0.5)
    return mp.sqrt(mp.pi / (2 * z)) * kind(order, z)


def _hankel_spherical(sign, n, z):
    if not _isint(n):
        return None
    order = n + mp.mpf(0.5)
    factor = mp.sqrt(mp.pi / (2 * z))
    return factor * (mp.besselj(order, z) + sign * 1j * mp.bessely(order, z))


def _meijerg(ap, bq, z, polar):
    # SymPy's meijerg._eval_evalf: the argument is z*exp_polar(polar). A
    # branch beyond (-pi, pi] is carried by mpmath's r argument, which
    # evaluates G(z'**(1/r)) with z' = z**(1/n)*exp(I*branch/n).
    if not polar:
        return mp.meijerg(ap, bq, z)
    branch = polar / 1j
    ratio = abs(branch / mp.pi)
    nearest = mp.nint(ratio)
    # An exact multiple of pi keeps SymPy's exact ceiling.
    if abs(ratio - nearest) <= ratio * mp.eps * 16:
        ratio = nearest
    n = int(mp.ceil(ratio)) + 1
    # mpmath forms the hypergeometric argument as z'**n from the rotated
    # root z'; on a sheet boundary (an odd multiple of pi) the rounding noise
    # of that power picks the side of the inner cut at random. The exact
    # power z*exp(I*branch) is passed instead, (-1)**k for branch = k*pi.
    if ratio == nearest:
        power = z if int(nearest) % 2 == 0 else -z
    else:
        power = z * mp.exp(1j * branch)
    z = z ** (mp.mpf(1) / n) * mp.exp(1j * branch / n)
    return _meijerg_power(ap, bq, z, mp.mpf(1) / n, power)


def _meijerg_power(a_s, b_s, z, r, power):
    # mpmath's meijerg(a_s, b_s, z, r) with z**(1/r) given exactly as power.
    an, ap = a_s
    bm, bq = b_s
    n = len(an)
    p = n + len(ap)
    m = len(bm)
    q = m + len(bq)
    a = [mp.convert(item) for item in an + ap]
    b = [mp.convert(item) for item in bm + bq]
    z = mp.convert(z)
    if p < q:
        series = 1
    elif p > q:
        series = 2
    else:
        series = 2 if m + n == p and abs(z) > 1 else 1
    if series == 1:
        def h(*args):
            a = args[:p]
            b = args[p:]
            terms = []
            for k in range(m):
                gn = [b[j] - b[k] for j in range(m) if j != k] + [1 - a[j] + b[k] for j in range(n)]
                gd = [a[j] - b[k] for j in range(n, p)] + [1 - b[j] + b[k] for j in range(m, q)]
                hn = [1 - a[j] + b[k] for j in range(p)]
                hd = [1 - b[j] + b[k] for j in range(q) if j != k]
                hz = (-mp.one) ** (p - m - n) * power
                terms.append(([z], [b[k] / r], gn, gd, hn, hd, hz))
            return terms
    else:
        def h(*args):
            a = args[:p]
            b = args[p:]
            terms = []
            for k in range(n):
                gn = [a[k] - a[j] for j in range(n) if j != k] + [1 - a[k] + b[j] for j in range(m)]
                gd = [a[k] - b[j] for j in range(m, q)] + [1 - a[k] + a[j] for j in range(n, p)]
                hn = [1 - a[k] + b[j] for j in range(q)]
                hd = [1 + a[j] - a[k] for j in range(p) if j != k]
                hz = (-mp.one) ** (q - m - n) / power
                terms.append(([z], [(a[k] - 1) / mp.convert(r)], gn, gd, hn, hd, hz))
            return terms
    return mp.hypercomb(h, a + b)


# name -> (function, extra argument precision, working precision offset)
_GENERIC = {
    "besselj": mp.besselj, "bessely": mp.bessely, "besseli": mp.besseli,
    "besselk": mp.besselk, "hankel1": mp.hankel1, "hankel2": mp.hankel2,
    "airyai": mp.airyai, "airybi": mp.airybi,
    "gamma": mp.gamma, "loggamma": mp.loggamma, "erf": mp.erf, "erfc": mp.erfc,
    "erfi": mp.erfi, "erfinv": mp.erfinv, "fresnels": mp.fresnels,
    "fresnelc": mp.fresnelc, "Ei": mp.ei, "li": mp.li, "Si": mp.si, "Ci": mp.ci,
    "Shi": mp.shi, "Chi": mp.chi, "expint": mp.expint,
    "zeta": mp.zeta, "dirichlet_eta": mp.altzeta, "polylog": mp.polylog,
    "lerchphi": mp.lerchphi, "beta": mp.beta,
    "elliptic_k": mp.ellipk, "elliptic_e": mp.ellipe, "elliptic_f": mp.ellipf,
    "elliptic_pi": _ellippi, "LambertW": mp.lambertw,
    "hyper": mp.hyper, "appellf1": mp.appellf1,
    "legendre": mp.legendre, "assoc_legendre": mp.legenp,
    "chebyshevt": mp.chebyt, "chebyshevu": mp.chebyu, "hermite": mp.hermite,
    "jacobi": mp.jacobi, "gegenbauer": mp.gegenbauer,
    # A non-integer degree evaluates the analytic continuation, as for
    # jacobi and gegenbauer (SymPy truncates the degree at construction).
    "assoc_laguerre": mp.laguerre,
    "factorial": mp.factorial, "binomial": mp.binomial,
    "RisingFactorial": mp.rf, "FallingFactorial": mp.ff,
}

_SPECIAL = {
    "lowergamma": (lambda a, z: mp.gammainc(a, 0, z), 0, 0),
    "uppergamma": (lambda a, z: mp.gammainc(a, z, mp.inf), 0, 0),
    "polygamma": (_polygamma, 12, 12),
    "digamma": (lambda z: _polygamma(mp.mpf(0), z), 12, 12),
    "trigamma": (lambda z: _polygamma(mp.mpf(1), z), 12, 12),
    "airyaiprime": (lambda z: mp.airyai(z, derivative=1), 0, 0),
    "airybiprime": (lambda z: mp.airybi(z, derivative=1), 0, 0),
    "Li": (lambda z: mp.li(z, offset=True), 5, 0),
    "jn": (lambda n, z: _spherical(mp.besselj, n, z), 5, 0),
    "yn": (lambda n, z: _spherical(mp.bessely, n, z), 5, 0),
    "hn1": (lambda n, z: _hankel_spherical(1, n, z), 5, 0),
    "hn2": (lambda n, z: _hankel_spherical(-1, n, z), 5, 0),
    "Ynm": (mp.spherharm, 0, 0),
    "bernoulli": (_bernoulli, 0, 0),
    "harmonic": (_harmonic, 0, 0),
    "euler": (_euler, 0, 0),
    "andre": (_andre, 12, 12),
    "meijerg": (_meijerg, 0, 0),
}


# Elementary heads whose mpmath value can be infinite at a Float argument.
_ELEMENTARY = {"atanh": mp.atanh, "acoth": mp.acoth, "atan": mp.atan, "acot": mp.acot}


def _code(raw):
    sign, man, exp, bc = raw
    if man:
        return 4
    if bc == 0 and exp == 0:
        return 0
    return {mpmath.libmp.finf: 1, mpmath.libmp.fninf: 2}.get(raw, 3)


def _nonfinite(name, args, prec):
    # SymPy's Function._eval_evalf converts mpmath's value with
    # Expr._from_mpmath: an infinite or nan part becomes oo, -oo or nan.
    # The parts are returned as codes (0 zero, 1 +inf, 2 -inf, 3 nan,
    # 4 finite nonzero) in the exponent fields of zero mantissas.
    if name in _SPECIAL:
        function, extra, offset = _SPECIAL[name]
    else:
        function = _GENERIC.get(name) or _ELEMENTARY.get(name)
        extra, offset = 5, 0
    if function is None:
        return None
    try:
        with workprec(prec + extra):
            values = [_convert(arg) for arg in args]
        with workprec(prec + offset):
            value = function(*values)
    except Exception:
        return None
    if isinstance(value, mpmath.mpf):
        codes = (_code(value._mpf_), 0)
    elif isinstance(value, mpmath.mpc):
        codes = tuple(_code(part) for part in value._mpc_)
    else:
        return None
    if not any(code in (1, 2, 3) for code in codes):
        return None
    return (0, codes[0], 0, codes[1])


def evaluate(name, args, prec):
    """Evaluate head ``name`` at ``prec`` bits, or return ``None``."""
    if name.startswith("!"):
        return _nonfinite(name[1:], args, prec)
    if name in _SPECIAL:
        function, extra, offset = _SPECIAL[name]
    elif name in _GENERIC:
        function, extra, offset = _GENERIC[name], 5, 0
    else:
        return None
    try:
        with workprec(prec + extra):
            values = [_convert(arg) for arg in args]
    except (ArithmeticError, ValueError, TypeError):
        return None
    # As in SymPy, an exception raised by mpmath itself propagates.
    with workprec(prec + offset):
        value = function(*values)
        if value is None:
            return None
        value = +value
    if offset:
        with workprec(prec):
            value = +value
    return _result(value)
