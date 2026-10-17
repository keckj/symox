"""Integral transforms with SymPy's API (``sympy.integrals.transforms`` and
``sympy.integrals.laplace``).

The transforms are computed natively by ``symox::transforms``; this module
only implements SymPy's signatures, hint handling (``noconds``,
``simplify``, ``needeval``, ``legacy_matrix``), the held transform classes
and their properties.
"""

import warnings

from . import _core

__all__ = [
    "IntegralTransform", "IntegralTransformError",
    "LaplaceTransform", "InverseLaplaceTransform",
    "FourierTransform", "InverseFourierTransform",
    "SineTransform", "InverseSineTransform",
    "CosineTransform", "InverseCosineTransform",
    "HankelTransform", "InverseHankelTransform",
    "MellinTransform", "InverseMellinTransform",
    "laplace_transform", "inverse_laplace_transform",
    "fourier_transform", "inverse_fourier_transform",
    "sine_transform", "inverse_sine_transform",
    "cosine_transform", "inverse_cosine_transform",
    "hankel_transform", "inverse_hankel_transform",
    "mellin_transform", "inverse_mellin_transform",
    "laplace_correspondence", "laplace_initial_conds",
]

_NONE = "_None"


class IntegralTransformError(NotImplementedError):
    """Exception raised when a transform cannot be computed (with the hint
    ``needeval=True``, or for Hankel transforms)."""

    def __init__(self, transform, function, msg):
        super().__init__(
            "%s Transform could not be computed: %s." % (transform, msg))
        self.function = function


def _sym():
    import symox
    return symox


def _is_matrix(f):
    return isinstance(f, _core.Matrix) or getattr(f, "_symox_sparse", False)


def _none_sentinel():
    return _sym().Symbol(_NONE)


def _from_sentinel(value):
    if isinstance(value, _core.Expr) and value.is_Symbol and value.name == _NONE:
        return None
    return value


# ---------------------------------------------------------------------------
# held transform classes
# ---------------------------------------------------------------------------

class _TransformMeta(type):
    def __instancecheck__(cls, value):
        if not isinstance(value, _core.Expr):
            return False
        name = _head_name(value)
        if name is None:
            return False
        if cls is IntegralTransform:
            return True
        return name == cls.__name__ or name in cls._subclasses

    def __subclasscheck__(cls, sub):
        return type.__subclasscheck__(cls, sub)


_HEADS = {
    "LaplaceTransform": ("Laplace", 3),
    "InverseLaplaceTransform": ("Inverse Laplace", 4),
    "FourierTransform": ("Fourier", 3),
    "InverseFourierTransform": ("Inverse Fourier", 3),
    "SineTransform": ("Sine", 3),
    "InverseSineTransform": ("Inverse Sine", 3),
    "CosineTransform": ("Cosine", 3),
    "InverseCosineTransform": ("Inverse Cosine", 3),
    "HankelTransform": ("Hankel", 4),
    "InverseHankelTransform": ("Inverse Hankel", 4),
    "MellinTransform": ("Mellin", 3),
    "InverseMellinTransform": ("Inverse Mellin", 5),
}


def _head_name(value):
    try:
        name = value.func.__name__
    except AttributeError:
        name = str(getattr(value, "func", ""))
    return name if name in _HEADS else None


class IntegralTransform(metaclass=_TransformMeta):
    """Base class of the unevaluated integral transforms."""
    _subclasses = ()
    _name = None

    def __new__(cls, *args, **options):
        if cls is IntegralTransform:
            raise TypeError("IntegralTransform is abstract")
        return _core._transform_held(cls.__name__, list(args))


class FourierTypeTransform(IntegralTransform):
    """Base class of the Fourier transforms."""
    _subclasses = ("FourierTransform", "InverseFourierTransform")


class SineCosineTypeTransform(IntegralTransform):
    """Base class of the sine and cosine transforms."""
    _subclasses = ("SineTransform", "InverseSineTransform",
                   "CosineTransform", "InverseCosineTransform")


class HankelTypeTransform(IntegralTransform):
    """Base class of the Hankel transforms."""
    _subclasses = ("HankelTransform", "InverseHankelTransform")


class LaplaceTransform(IntegralTransform):
    """Unevaluated Laplace transform ``LaplaceTransform(f, t, s)``."""
    _name = "Laplace"


class InverseLaplaceTransform(IntegralTransform):
    """Unevaluated inverse Laplace transform
    ``InverseLaplaceTransform(F, s, t, plane)``."""
    _name = "Inverse Laplace"

    def __new__(cls, F, s, x, plane, **opts):
        if plane is None:
            plane = _none_sentinel()
        return _core._transform_held(cls.__name__, [F, s, x, plane])


class FourierTransform(FourierTypeTransform):
    """Unevaluated Fourier transform ``FourierTransform(f, x, k)``."""
    _name = "Fourier"
    _subclasses = ()


class InverseFourierTransform(FourierTypeTransform):
    """Unevaluated inverse Fourier transform."""
    _name = "Inverse Fourier"
    _subclasses = ()


class SineTransform(SineCosineTypeTransform):
    """Unevaluated sine transform."""
    _name = "Sine"
    _subclasses = ()


class InverseSineTransform(SineCosineTypeTransform):
    """Unevaluated inverse sine transform."""
    _name = "Inverse Sine"
    _subclasses = ()


class CosineTransform(SineCosineTypeTransform):
    """Unevaluated cosine transform."""
    _name = "Cosine"
    _subclasses = ()


class InverseCosineTransform(SineCosineTypeTransform):
    """Unevaluated inverse cosine transform."""
    _name = "Inverse Cosine"
    _subclasses = ()


class HankelTransform(HankelTypeTransform):
    """Unevaluated Hankel transform ``HankelTransform(f, r, k, nu)``."""
    _name = "Hankel"
    _subclasses = ()


class InverseHankelTransform(HankelTypeTransform):
    """Unevaluated inverse Hankel transform."""
    _name = "Inverse Hankel"
    _subclasses = ()


class MellinTransform(IntegralTransform):
    """Unevaluated Mellin transform ``MellinTransform(f, x, s)``."""
    _name = "Mellin"


class InverseMellinTransform(IntegralTransform):
    """Unevaluated inverse Mellin transform
    ``InverseMellinTransform(F, s, x, a, b)``."""
    _name = "Inverse Mellin"

    def __new__(cls, F, s, x, a, b, **opts):
        a = _none_sentinel() if a is None else a
        b = _none_sentinel() if b is None else b
        return _core._transform_held(cls.__name__, [F, s, x, a, b])


# ---------------------------------------------------------------------------
# doit and properties (called by the native Expr)
# ---------------------------------------------------------------------------

# (noconds, simplify) defaults of each class's doit
_DOIT_DEFAULTS = {
    "LaplaceTransform": (True, False),
    "InverseLaplaceTransform": (True, False),
    "MellinTransform": (False, True),
    "InverseMellinTransform": (True, True),
}


def _shape(name, items):
    """Arrange the native result list like SymPy's return value."""
    if len(items) == 1:
        return items[0]
    if name == "MellinTransform":
        return items[0], (items[1], items[2]), items[3]
    return tuple(items)


def _doit(expr, deep=True, **hints):
    name = _head_name(expr)
    noconds_default, simplify_default = _DOIT_DEFAULTS.get(name, (True, True))
    noconds = hints.pop("noconds", noconds_default)
    simplify = hints.pop("simplify", simplify_default)
    needeval = hints.pop("needeval", False)
    hints.pop("legacy_matrix", None)
    args = list(expr.args)
    if deep:
        args = [a.doit() if isinstance(a, _core.Expr) else a for a in args]
    try:
        items = _core._transform_doit(name, args, bool(simplify), bool(needeval))
    except NotImplementedError as err:
        msg = str(err).split(": ", 1)[-1]
        raise IntegralTransformError(_HEADS[name][0], args[0], msg) from None
    if noconds:
        return items[0]
    return _shape(name, items)


def _as_integral(expr):
    s = _sym()
    name = _head_name(expr)
    a = expr.args
    I, oo, pi = s.I, s.oo, s.pi
    if name == "LaplaceTransform":
        f, t, p = a
        return s.Integral(f*s.exp(-p*t), (t, 0, oo))
    if name == "InverseLaplaceTransform":
        F, p, t = a[:3]
        c = s.Dummy("c")
        return s.Integral(s.exp(p*t)*F, (p, c - I*oo, c + I*oo))/(2*pi*I)
    if name == "MellinTransform":
        f, x, p = a
        return s.Integral(f*x**(p - 1), (x, 0, oo))
    if name == "InverseMellinTransform":
        F, p, x = a[:3]
        c = s.Dummy("c")
        return s.Integral(F*x**(-p), (p, c - I*oo, c + I*oo))/(2*pi*I)
    if name in ("FourierTransform", "InverseFourierTransform"):
        f, x, k = a
        b = -2*pi if name == "FourierTransform" else 2*pi
        return s.Integral(f*s.exp(b*I*x*k), (x, -oo, oo))
    if name in SineCosineTypeTransform._subclasses:
        f, x, k = a
        K = s.sin if "Sine" in name else s.cos
        return s.Integral(s.sqrt(2)/s.sqrt(pi)*f*K(x*k), (x, 0, oo))
    if name in HankelTypeTransform._subclasses:
        f, r, k, nu = a
        return s.Integral(f*s.besselj(nu, k*r)*r, (r, 0, oo))
    raise AttributeError("as_integral")


def _attribute(expr, attr):
    name = _head_name(expr)
    if attr == "as_integral":
        return _as_integral(expr)
    if attr == "fundamental_plane" and name == "InverseLaplaceTransform":
        return _from_sentinel(expr.args[3])
    if attr == "fundamental_strip" and name == "InverseMellinTransform":
        return _from_sentinel(expr.args[3]), _from_sentinel(expr.args[4])
    raise AttributeError("%s has no attribute %s" % (name, attr))


# ---------------------------------------------------------------------------
# public functions
# ---------------------------------------------------------------------------

def _apply_matrix(f, fn):
    return f.applyfunc(fn)


def laplace_transform(f, t, s, legacy_matrix=True, **hints):
    """Laplace transform ``F(s) = int_0^oo f(t) exp(-s*t) dt``.

    Returns ``(F, a, cond)`` where ``Re(s) > a`` is the half plane of
    convergence, or ``F`` alone with ``noconds=True``.
    """
    noconds = hints.get("noconds", False)
    simplify = hints.get("simplify", False)
    if _is_matrix(f) and hasattr(f, "applyfunc"):
        conds = not noconds
        if conds and legacy_matrix:
            warnings.warn(
                "Calling laplace_transform() on a Matrix with noconds=False "
                "(the default) is deprecated. Either noconds=True or use "
                "legacy_matrix=False to get the new behavior.",
                DeprecationWarning, stacklevel=2)
            return f.applyfunc(
                lambda fij: _sym().Tuple(*laplace_transform(fij, t, s, **hints)))
        elements = [laplace_transform(fij, t, s, **hints) for fij in f]
        if conds:
            s_ = _sym()
            els, avals, cnds = zip(*elements)
            return (s_.Matrix(f.shape[0], f.shape[1], list(els)),
                    s_.Max(*avals), s_.And(*cnds))
        return _sym().Matrix(f.shape[0], f.shape[1], elements)
    r = _core._transform_doit("LaplaceTransform", [f, t, s], bool(simplify), False)
    if noconds:
        return r[0]
    return tuple(r)


def inverse_laplace_transform(F, s, t, plane=None, **hints):
    """Inverse Laplace transform of ``F(s)``; ``(f, cond)`` with
    ``noconds=False``."""
    noconds = hints.get("noconds", True)
    simplify = hints.get("simplify", False)
    if _is_matrix(F) and hasattr(F, "applyfunc"):
        return F.applyfunc(
            lambda Fij: inverse_laplace_transform(Fij, s, t, plane, **hints))
    if plane is None:
        plane = _none_sentinel()
    r = _core._transform_doit("InverseLaplaceTransform", [F, s, t, plane], bool(simplify), False)
    if noconds:
        return r[0]
    return r[0], r[1]


def _generic(name, args, hints):
    return _doit(_core._transform_held(name, list(args)), deep=False, **hints)


def mellin_transform(f, x, s, **hints):
    """Mellin transform ``F(s) = int_0^oo x**(s-1) f(x) dx``: returns
    ``(F, (a, b), cond)`` where ``a < Re(s) < b`` is the fundamental strip."""
    return _generic("MellinTransform", (f, x, s), hints)


def inverse_mellin_transform(F, s, x, strip, **hints):
    """Inverse Mellin transform of ``F(s)`` on the strip ``(a, b)``."""
    a, b = strip
    return _doit(InverseMellinTransform(F, s, x, a, b), deep=False, **hints)


def fourier_transform(f, x, k, **hints):
    """Fourier transform ``F(k) = int f(x) exp(-2*pi*I*x*k) dx``."""
    return _generic("FourierTransform", (f, x, k), hints)


def inverse_fourier_transform(F, k, x, **hints):
    """Inverse Fourier transform ``f(x) = int F(k) exp(2*pi*I*x*k) dk``."""
    return _generic("InverseFourierTransform", (F, k, x), hints)


def sine_transform(f, x, k, **hints):
    """Sine transform ``sqrt(2/pi) int_0^oo f(x) sin(k*x) dx``."""
    return _generic("SineTransform", (f, x, k), hints)


def inverse_sine_transform(F, k, x, **hints):
    """Inverse sine transform ``sqrt(2/pi) int_0^oo F(k) sin(k*x) dk``."""
    return _generic("InverseSineTransform", (F, k, x), hints)


def cosine_transform(f, x, k, **hints):
    """Cosine transform ``sqrt(2/pi) int_0^oo f(x) cos(k*x) dx``."""
    return _generic("CosineTransform", (f, x, k), hints)


def inverse_cosine_transform(F, k, x, **hints):
    """Inverse cosine transform ``sqrt(2/pi) int_0^oo F(k) cos(k*x) dk``."""
    return _generic("InverseCosineTransform", (F, k, x), hints)


def hankel_transform(f, r, k, nu, **hints):
    """Hankel transform ``F_nu(k) = int_0^oo f(r) J_nu(k*r) r dr``."""
    return _generic("HankelTransform", (f, r, k, nu), hints)


def inverse_hankel_transform(F, k, r, nu, **hints):
    """Inverse Hankel transform ``f(r) = int_0^oo F(k) J_nu(k*r) k dk``."""
    return _generic("InverseHankelTransform", (F, k, r, nu), hints)


# ---------------------------------------------------------------------------
# Laplace helpers
# ---------------------------------------------------------------------------

def _rebuild(f, fn):
    args = f.args
    if not args:
        return f
    new = [fn(a) if isinstance(a, _core.Expr) else a for a in args]
    if all(n is o for n, o in zip(new, args)):
        return f
    return f.func(*new)


def laplace_correspondence(f, fdict, /):
    """Replace ``LaplaceTransform(y(t), t, s)`` by ``Y(s)`` and
    ``InverseLaplaceTransform(Y(s), s, t, p)`` by ``y(t)`` for each pair
    ``y: Y`` of ``fdict``."""
    if not isinstance(f, _core.Expr):
        return f
    if not (f.has(LaplaceTransform) if False else _has_laplace(f)):
        return f
    name = _head_name(f)
    for y, Y in fdict.items():
        if name == "LaplaceTransform":
            g, t, s = f.args
            if _applied(g, y) and g.args[0] == t:
                return Y(s)
        if name == "InverseLaplaceTransform":
            G, s, t = f.args[:3]
            if _applied(G, Y) and G.args[0] == s:
                return y(t)
    return _rebuild(f, lambda a: laplace_correspondence(a, fdict))


def _applied(g, y):
    try:
        return len(g.args) == 1 and g.func == y
    except Exception:
        return False


def _has_laplace(f):
    found = []

    def walk(e):
        if found:
            return
        if _head_name(e) in ("LaplaceTransform", "InverseLaplaceTransform"):
            found.append(e)
            return
        for a in e.args:
            if isinstance(a, _core.Expr):
                walk(a)
    walk(f)
    return bool(found)


def laplace_initial_conds(f, t, fdict, /):
    """Substitute the initial conditions ``fdict = {y: [y(0), y'(0), ...]}``
    into the result of a Laplace transform."""
    s = _sym()
    for y, ic in fdict.items():
        for k in range(len(ic)):
            if k == 0:
                old = y(0)
            elif k == 1:
                old = s.Subs(s.Derivative(y(t), t), t, 0)
            else:
                old = s.Subs(s.Derivative(y(t), (t, k)), t, 0)
            f = f.xreplace({old: ic[k]})
    return f
