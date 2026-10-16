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
    Wild,
    __version__,
    cbrt,
    diff,
    expand,
    integrate,
    lambdify,
    log,
    radsimp,
    root,
    simplify,
    sqrt,
    srepr,
    summation,
    bernoulli,
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
# number theory
# ---------------------------------------------------------------------------

from ._core import (  # noqa: E402
    binomial, crt, divisor_count, divisor_sigma, divisors, factorial, factorint,
    integer_nthroot, isprime, jacobi_symbol, legendre_symbol, mobius, mod_inverse,
    n_order, nextprime, perfect_power, prevprime, prime, primepi, primerange,
    primitive_root, totient,
)

from ._core import (  # noqa: E402
    expand_log, expand_mul, expand_power_base, expand_power_exp, expand_trig,
    logcombine, nsimplify, powsimp, trigsimp,
)
from ._core import checkodesol, dsolve, nsolve, solve  # noqa: E402

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

def _expand_range(name):
    """``x:3`` -> x0 x1 x2, ``x1:3`` -> x1 x2, ``a:d`` -> a b c d."""
    if ":" not in name:
        return [name]
    left, right = name.split(":", 1)
    if right.isdigit() or (right == "" and left[-1:].isdigit()):
        if right == "":
            raise ValueError("missing end of range")
        m = _re.match(r"^(.*?)(\d*)$", left)
        prefix, start = m.group(1), m.group(2)
        return [f"{prefix}{i}" for i in range(int(start) if start else 0, int(right))]
    letters = _string.ascii_letters
    if left and left[-1] in letters:
        prefix, start = left[:-1], left[-1]
    else:
        prefix, start = left, "a"
    if right and right not in letters:
        raise ValueError(f"invalid range in symbol name {name!r}")
    lo = letters.index(start)
    hi = letters.index(right) + 1 if right else len(letters)
    return [f"{prefix}{c}" for c in letters[lo:hi]]


def symbols(names, *, seq=False, cls=None, **assumptions):
    if cls is Wild:
        syms = symbols(names, seq=True)
        out = tuple(Wild(s.name, **assumptions) for s in syms)
        return out if (seq or len(out) != 1) else out[0]
    return _symbols(names, seq=seq, cls=cls, **assumptions)


def _symbols(names, *, seq=False, cls=None, **assumptions):
    """Create symbols like ``sympy.symbols`` (``"x y"``, ``"x,y"``, ``"x:3"``)."""
    if not isinstance(names, str):
        return type(names)(symbols(n, cls=cls, **assumptions) for n in names)
    parts = [p for p in _re.split(r"[\s,]+", names.strip()) if p]
    out = []
    for p in parts:
        for n in _expand_range(p):
            out.append(Function(n) if cls is Function else Symbol(n, **assumptions))
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


# --- poly ---
from ._core import (  # noqa: E402
    Poly,
    apart,
    cancel,
    cofactors,
    collect,
    compose,
    decompose,
    degree,
    discriminant,
    div,
    exquo,
    factor,
    factor_list,
    gcd,
    gcdex,
    half_gcdex,
    invert,
    lcm,
    nroots,
    pdiv,
    pexquo,
    pquo,
    prem,
    quo,
    real_roots,
    reduced,
    rem,
    resultant,
    roots,
    sqf,
    sqf_list,
    sqf_part,
    terms_gcd,
)

together = _core._poly_together


class PolynomialError(ValueError):
    """Raised for invalid polynomial operations."""


def gcd_list(seq, *gens, **args):
    return gcd(list(seq), *gens, **args)


def lcm_list(seq, *gens, **args):
    return lcm(list(seq), *gens, **args)


def poly(expr, *gens, **args):
    return Poly(expr, *gens, **args)


def LC(f, *gens, **args):
    """Leading coefficient of ``f``."""
    return Poly(f, *gens, **args).LC()


def LM(f, *gens, **args):
    """Leading monomial of ``f`` (as an expression)."""
    p = Poly(f, *gens, **args)
    m = 1
    for g, k in zip(p.gens, p.LM()):
        m = m * g**k
    return sympify(m)


def LT(f, *gens, **args):
    """Leading term of ``f``."""
    p = Poly(f, *gens, **args)
    return p.LC() * LM(f, *gens, **args)


def primitive(f, *gens, **args):
    """``(content, primitive part)`` of a polynomial."""
    c, p = Poly(f, *gens, **args).primitive()
    return c, p.as_expr()


def _install_poly_methods():
    methods = {
        "factor": lambda self, *gens, **kw: factor(self, *gens, **kw),
        "cancel": lambda self, *gens, **kw: cancel(self, *gens, **kw),
        "apart": lambda self, x=None, **kw: apart(self, x, **kw),
        "collect": lambda self, syms, **kw: collect(self, syms, **kw),
        "as_poly": lambda self, *gens, **kw: Poly(self, *gens, **kw),
        "together": lambda self, *a, **kw: together(self, *a, **kw),
    }
    for name, fn in methods.items():
        try:
            setattr(Expr, name, fn)
        except (TypeError, AttributeError):  # immutable extension type
            pass
    try:
        Poly.get_domain = lambda self: self.domain

        def _get_modulus(self):
            if self.modulus is None:
                raise PolynomialError("not a polynomial over a Galois field")
            return self.modulus

        Poly.get_modulus = _get_modulus
    except (TypeError, AttributeError):
        pass


_install_poly_methods()
# --- end poly ---
# --- series ---
# Series expansions, the Order term and limits (crates/symox-py/src/series.rs).


class PoleError(ValueError):
    """Raised when an expansion point is a singularity of the expression."""


def _series_call(expr, x, x0, n, dir):
    try:
        return _core._series(expr, x, x0, n, dir)
    except ValueError as e:
        msg = str(e)
        if msg.startswith("PoleError: "):
            raise PoleError(msg[len("PoleError: "):]) from None
        raise


def series(expr, x=None, x0=0, n=6, dir="+"):
    """Series expansion of ``expr`` around ``x = x0`` up to ``O((x - x0)**n)``."""
    return _series_call(expr, x, x0, n, dir)


def Order(expr, *args, **kwargs):
    """Landau order term ``O(expr, (x, x0), ...)``."""
    return _core._order(expr, *args)


O = Order


def _expr_series(self, x=None, x0=0, n=6, dir="+", logx=None, cdir=0):
    return _series_call(self, x, x0, n, dir)


def _expr_nseries(self, x=None, x0=0, n=6, dir="+", logx=None, cdir=0):
    return _series_call(self, x, x0, n, dir)


def _expr_taylor_term(self, n, x, *previous_terms):
    return _core._taylor_term(self, n, x)


def _expr_getO(self):
    return _core._get_o(self)


def _expr_getn(self):
    n = _core._getn(self)
    return n


def _expr_lseries(self, x=None, x0=0, dir="+", logx=None, cdir=0):
    """Lazy iterator over the terms of the expansion (in increasing order)."""
    n = 6
    done = 0
    while True:
        s = _series_call(self, x, x0, n, dir)
        terms = s.removeO()
        terms = list(terms.args) if terms.is_Add else ([terms] if terms != 0 else [])
        terms = _sort_series_terms(terms, s)
        for t in terms[done:]:
            yield t
        done = len(terms)
        if _core._get_o(s) is None:
            return
        n *= 2


def _sort_series_terms(terms, s):
    order = [str(t) for t in (s.as_ordered_terms() if s.is_Add else [s])]
    return sorted(terms, key=lambda t: order.index(str(t)) if str(t) in order else len(order))


def _expr_leadterm(self, x, x0=0, dir="+"):
    return _core._leadterm(self, x, x0, dir)


def _expr_as_leading_term(self, *symbols, **kwargs):
    x = symbols[0] if symbols else None
    c, e = _core._leadterm(self, x)
    return c * x**e


Expr.series = _expr_series
Expr.nseries = _expr_nseries
Expr.taylor_term = _expr_taylor_term
Expr.removeO = lambda self: _core._remove_o(self)
Expr.getO = _expr_getO
Expr.getn = _expr_getn
Expr.lseries = _expr_lseries
Expr.leadterm = _expr_leadterm
Expr.as_leading_term = _expr_as_leading_term
Expr.is_Order = property(lambda self: _core._is_order(self))
Expr.contains = lambda self, t: _core._order_contains(self, t)


def _order_part(i):
    def get(self):
        if not _core._is_order(self):
            raise AttributeError("only Order terms have expr/variables/point")
        a = self.args
        if i == 0:
            return a[0]
        return tuple(a[1::2]) if i == 1 else tuple(a[2::2])
    return property(get)


Expr.expr = _order_part(0)
Expr.variables = _order_part(1)
Expr.point = _order_part(2)

class Limit:
    """Unevaluated limit ``Limit(e, z, z0, dir)``; ``doit()`` evaluates it."""

    def __init__(self, e, z, z0, dir="+"):
        self.expr = sympify(e)
        self.var = sympify(z)
        self.point = sympify(z0)
        if str(self.point) == "oo":
            dir = "-"
        elif str(self.point) == "-oo":
            dir = "+"
        self.dir = dir
        self.args = (self.expr, self.var, self.point, dir)

    def doit(self, **hints):
        r = _core._limit(self.expr, self.var, self.point, self.dir)
        return self if r is None else r

    def __str__(self):
        return f"Limit({self.expr}, {self.var}, {self.point}, dir='{self.dir}')"

    __repr__ = __str__

    def __eq__(self, other):
        return isinstance(other, Limit) and self.args == other.args

    def __hash__(self):
        return hash(("Limit",) + tuple(str(a) for a in self.args))


def limit(e, z, z0, dir="+"):
    """Limit of ``e`` as ``z -> z0`` (``dir`` is ``'+'``, ``'-'`` or ``'+-'``).

    Returns an unevaluated :class:`Limit` when the limit cannot be computed.
    """
    r = _core._limit(e, z, z0, dir)
    return Limit(e, z, z0, dir) if r is None else r


Expr.limit = lambda self, x, xlim, dir="+": limit(self, x, xlim, dir)
# --- end series ---

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

# SymPy's constructor classes and printer extension points. Keep this last:
# constructors delegate to the native API defined above.
from ._compat import (  # noqa: E402
    Add, Mul, Pow, Symbol, Number, Integer, Rational, Float, UnevaluatedExpr,
    S, count_ops, default_sort_key,
)
from .codegen.ast import Indexed, IndexedBase  # noqa: E402
from .printing.str import StrPrinter  # noqa: E402
from .printing.pretty import PrettyPrinter  # noqa: E402

__all__ = [n for n in dir() if not n.startswith("_")]
