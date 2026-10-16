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


# --- printing ---
# LaTeX / pretty printers, common subexpression elimination and code
# generation (C, C++, CUDA, Python, NumPy, Rust, Fortran).
import builtins as _builtins  # noqa: E402

from ._core import (  # noqa: E402
    ccode,
    ccode_function,
    cuda_device_function,
    cuda_kernel,
    cudacode,
    cxxcode,
    fcode,
    latex,
    numpycode,
    pretty,
    pycode,
    rust_code,
)


def pprint(expr, **settings):
    """Print the 2D unicode form of ``expr`` (SymPy's ``pprint``)."""
    print(pretty(expr, **settings))


pretty_print = pprint


def print_latex(expr, **settings):
    print(latex(expr, **settings))


def print_ccode(expr, **settings):
    print(ccode(expr, **settings))


def numbered_symbols(prefix="x", start=0, exclude=()):
    """Infinite generator of symbols ``prefix0, prefix1, ...``."""
    excluded = {str(e) for e in exclude}
    i = start
    while True:
        name = f"{prefix}{i}"
        i += 1
        if name not in excluded:
            yield Symbol(name)


def cse(exprs, symbols=None, optimizations=None, postprocess=None,
        order="canonical", ignore=(), list=True):
    """Common subexpression elimination (SymPy's ``cse``).

    Returns ``(replacements, reduced_exprs)`` where ``replacements`` is a
    list of ``(Symbol, expr)`` pairs in evaluation order.
    """
    single = isinstance(exprs, (Expr, int, float, str))
    if isinstance(exprs, dict):
        keys = _builtins.list(exprs)
        reps, red = cse([exprs[k] for k in keys], symbols, optimizations, postprocess, order, ignore)
        return reps, dict(zip(keys, red))
    seq = [exprs] if single else exprs
    opt = optimizations not in (False, "none")
    reps, red = _core._cse(seq, opt, ignore or None)
    if symbols is not None:
        it = iter(symbols)
        mapping = {}
        new_reps = []
        for s, v in reps:
            v = v.xreplace(mapping) if mapping else v
            t = next(it)
            mapping[s] = t
            new_reps.append((t, v))
        reps = new_reps
        red = [r.xreplace(mapping) for r in red] if mapping else red
    if postprocess is not None:
        return postprocess(reps, red)
    if single and not list:
        return reps, red[0]
    if not single and not list:
        return reps, type(exprs)(red)
    return reps, red

# --- end printing ---

__all__ = [n for n in dir() if not n.startswith("_")]
