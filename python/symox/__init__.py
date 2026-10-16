"""symox: fast symbolic mathematics in Rust with a SymPy compatible API.

>>> from symox import symbols, sin, expand, diff
>>> x, y = symbols("x y")
>>> expand((x + y)**2)
x**2 + 2*x*y + y**2
>>> diff(sin(x)*x, x)
x*cos(x) + sin(x)
"""

from __future__ import annotations

import re as _re
import string as _string

from . import _core
from ._core import (
    CompiledFunction,
    Derivative,
    Expr,
    Float,
    Integer,
    N,
    Piecewise,
    Rational,
    Symbol,
    __version__,
    cbrt,
    diff,
    expand,
    lambdify,
    log,
    radsimp,
    root,
    simplify,
    sqrt,
    srepr,
    sympify,
    together,
)

Basic = Expr
S = sympify

# ---------------------------------------------------------------------------
# constants
# ---------------------------------------------------------------------------

pi = _core.constant("pi")
E = _core.constant("E")
I = _core.constant("I")
oo = _core.constant("oo")
zoo = _core.constant("zoo")
nan = _core.constant("nan")
EulerGamma = _core.constant("EulerGamma")
true = _core.constant("true")
false = _core.constant("false")

# ---------------------------------------------------------------------------
# functions
# ---------------------------------------------------------------------------


def _make(name):
    def f(*args):
        return _core._apply(name, *args)

    f.__name__ = name
    f.__qualname__ = name
    f.__doc__ = f"symbolic {name}"
    return f


_FUNCS = [
    "exp", "sin", "cos", "tan", "cot", "sec", "csc", "asin", "acos", "atan",
    "acot", "atan2", "sinh", "cosh", "tanh", "coth", "asinh", "acosh", "atanh",
    "Abs", "sign", "floor", "ceiling", "gamma", "loggamma", "erf", "erfc",
    "re", "im", "conjugate", "arg", "Heaviside", "DiracDelta", "Max", "Min",
]
for _n in _FUNCS:
    globals()[_n] = _make(_n)
ln = log
Abs = globals()["Abs"]


class Function:
    """Undefined function: ``f = Function('f'); f(x)``."""

    def __init__(self, name: str, **_assumptions):
        self.name = name
        self.__name__ = name

    def __call__(self, *args):
        return _core._apply(self.name, *args)

    def __repr__(self):
        return self.name

    def __eq__(self, other):
        return isinstance(other, Function) and other.name == self.name

    def __hash__(self):
        return hash(("Function", self.name))


# ---------------------------------------------------------------------------
# relationals
# ---------------------------------------------------------------------------


def Eq(a, b=0):
    return _core.Rel("==", a, b)


def Ne(a, b):
    return _core.Rel("!=", a, b)


def Lt(a, b):
    return _core.Rel("<", a, b)


def Le(a, b):
    return _core.Rel("<=", a, b)


def Gt(a, b):
    return _core.Rel(">", a, b)


def Ge(a, b):
    return _core.Rel(">=", a, b)


# ---------------------------------------------------------------------------
# symbols()
# ---------------------------------------------------------------------------

_RANGE = _re.compile(r"^(.*?)(\d*|[a-zA-Z]?):(\d*|[a-zA-Z]?)$")


def _expand_range(name):
    if ":" not in name:
        return [name]
    m = _RANGE.match(name)
    if not m:
        raise ValueError(f"invalid range in symbol name {name!r}")
    prefix, a, b = m.groups()
    if (a == "" or a.isdigit()) and (b == "" or b.isdigit()):
        if b == "":
            raise ValueError("missing end of range")
        lo = int(a) if a else 0
        return [f"{prefix}{i}" for i in range(lo, int(b))]
    letters = _string.ascii_letters
    lo = letters.index(a) if a else 0
    hi = letters.index(b) + 1 if b else len(letters)
    return [f"{prefix}{c}" for c in letters[lo:hi]]


def symbols(names, *, seq=False, cls=None, **assumptions):
    """Create symbols like ``sympy.symbols`` (``"x y"``, ``"x,y"``, ``"x:3"``)."""
    if not isinstance(names, str):
        return type(names)(symbols(n, **assumptions) for n in names)
    parts = [p for p in _re.split(r"[\s,]+", names.strip()) if p]
    out = []
    for p in parts:
        for n in _expand_range(p):
            out.append(Function(n) if cls is Function else Symbol(n))
    if len(out) == 1 and not seq and "," not in names and ":" not in names:
        return out[0]
    return tuple(out)


def var(names, **assumptions):
    """Create symbols and inject them in the caller's global namespace."""
    import inspect

    syms = symbols(names, seq=True, **assumptions)
    frame = inspect.currentframe().f_back
    for s in syms:
        frame.f_globals[s.name] = s
    return syms[0] if len(syms) == 1 else syms


# ---------------------------------------------------------------------------
# SymPy interop
# ---------------------------------------------------------------------------


def from_sympy(e):
    """Convert a SymPy expression (or anything sympifiable) to symox."""
    return sympify(e)


def to_sympy(e):
    """Convert a symox expression to SymPy."""
    return sympify(e).to_sympy()


# --- matrices ---
# Symbolic matrices (Rust implementation in symox._core.Matrix).

from ._core import Matrix, diag, eye, hessian, jacobian, ones, wronskian, zeros  # noqa: E402

MutableDenseMatrix = Matrix
DenseMatrix = Matrix
ImmutableMatrix = Matrix


class MatrixError(Exception):
    """Base class for matrix errors (as in SymPy)."""


class ShapeError(ValueError, MatrixError):
    """Wrong matrix shape."""


class NonSquareMatrixError(ShapeError):
    """The matrix must be square."""


class NonInvertibleMatrixError(ValueError, MatrixError):
    """The matrix is singular."""


class PurePoly:
    """Minimal ``PurePoly`` returned by ``Matrix.charpoly``."""

    __slots__ = ("_coeffs", "_gen", "_domain", "_expr")

    @classmethod
    def _from_coeffs(cls, coeffs, gen, domain, expr):
        p = cls.__new__(cls)
        p._coeffs = list(coeffs)
        p._gen = gen
        p._domain = domain
        p._expr = expr
        return p

    def as_expr(self, *gens):
        if gens:
            return self._expr.subs(self._gen, gens[0])
        return self._expr

    def all_coeffs(self):
        return list(self._coeffs)

    def coeffs(self):
        return [c for c in self._coeffs if c != 0]

    def degree(self, *gens):
        return len(self._coeffs) - 1

    @property
    def gen(self):
        return self._gen

    @property
    def gens(self):
        return (self._gen,)

    @property
    def domain(self):
        return self._domain

    def get_domain(self):
        return self._domain

    def eval(self, x):
        return self._expr.subs(self._gen, x)

    __call__ = eval

    def __eq__(self, other):
        if isinstance(other, PurePoly):
            return self._coeffs == other._coeffs
        return self._expr == other

    def __hash__(self):
        return hash(tuple(self._coeffs))

    def __str__(self):
        g = str(self._gen)
        n = len(self._coeffs) - 1
        terms = []
        for i, c in enumerate(self._coeffs):
            e = n - i
            if c == 0:
                continue
            mon = "" if e == 0 else (g if e == 1 else "%s**%d" % (g, e))
            if c.is_Add:
                sc = "(" + str(c) + ")" if mon else str(c)
            else:
                if mon and c == 1:
                    terms += ["+", mon]
                    continue
                if mon and c == -1:
                    terms += ["-", mon]
                    continue
                sc = str(c)
            st = sc if not mon else sc + "*" + mon
            if st.startswith("-"):
                terms += ["-", st[1:]]
            else:
                terms += ["+", st]
        if not terms:
            terms = ["+", "0"]
        if terms[0] in ("-", "+"):
            mod = terms.pop(0)
            if mod == "-":
                terms[0] = "-" + terms[0]
        return "PurePoly(%s, %s, domain='%s')" % (" ".join(terms), g, self._domain)

    __repr__ = __str__


def det(M, **kwargs):
    """Determinant of a matrix."""
    return Matrix(M).det(**kwargs)


def trace(M):
    return Matrix(M).trace()


def transpose(M):
    return Matrix(M).T


# --- end matrices ---

__all__ = [n for n in dir() if not n.startswith("_")]
