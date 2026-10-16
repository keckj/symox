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
    crt, divisor_count, divisor_sigma, divisors, factorint,
    integer_nthroot, isprime, jacobi_symbol, legendre_symbol, mobius, mod_inverse,
    n_order, nextprime, perfect_power, prevprime, prime, primepi, primerange,
    primitive_root, totient,
)

from ._core import (  # noqa: E402
    expand_log, expand_mul, expand_power_base, expand_power_exp, expand_trig,
    logcombine, nsimplify, powsimp, trigsimp,
)
from ._core import checkodesol, dsolve, nsolve, solve  # noqa: E402
from ._core import NonlinearError  # noqa: E402
from ._linear import linear_eq_to_matrix, linsolve, solve_linear_system  # noqa: E402
from ._polysolve import solve_poly_system, nonlinsolve  # noqa: E402
from ._core import residue  # noqa: E402

# --- special ---
# Special functions (symox::special): every head is a real function node with
# SymPy's automatic evaluation, derivatives and numerics.
from ._core import (  # noqa: E402
    besselsimp, expand_func, hyperexpand, multinomial_coefficients, stirling,
)

_SPECIAL = [
    "factorial", "factorial2", "subfactorial", "binomial", "RisingFactorial",
    "FallingFactorial", "harmonic", "bernoulli", "euler", "catalan", "fibonacci",
    "lucas", "partition", "polygamma", "digamma", "trigamma", "beta",
    "lowergamma", "uppergamma", "zeta", "dirichlet_eta", "polylog", "lerchphi",
    "erfi", "erf2", "erfinv", "erfcinv", "fresnels", "fresnelc", "Ei", "li", "Li",
    "Si", "Ci", "Shi", "Chi", "expint", "E1", "besselj", "bessely", "besseli",
    "besselk", "hankel1", "hankel2", "jn", "yn", "airyai", "airybi",
    "airyaiprime", "airybiprime", "elliptic_k", "elliptic_e", "elliptic_f",
    "elliptic_pi", "LambertW", "legendre", "assoc_legendre", "chebyshevt",
    "chebyshevu", "hermite", "laguerre", "assoc_laguerre", "jacobi",
    "gegenbauer", "Ynm",
]
for _n in _SPECIAL:
    globals()[_n] = _make(_n)
rf = RisingFactorial  # noqa: F821
ff = FallingFactorial  # noqa: F821


def _tuple(items):
    return _core._apply("Tuple", *items)


def hyper(ap, bq, z):
    """Generalised hypergeometric function ``hyper(ap, bq, z)``."""
    return _core._apply("hyper", _tuple(ap), _tuple(bq), z)


def bell(n, k_sym=None, symbols=None):
    """Bell numbers, Bell polynomials and incomplete Bell polynomials."""
    if k_sym is None:
        return _core._apply("bell", n)
    if symbols is None:
        return _core._apply("bell", n, k_sym)
    return _core._apply("bell", n, k_sym, _tuple(symbols))


def special_nstr(e, dps=15):
    """Arbitrary precision value of ``e`` as a string (``None`` if not
    supported)."""
    return _core._nstr(e, dps)
# --- end special ---

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


# --- matrices ---
# Symbolic matrices (Rust implementation in symox._core.Matrix).

from ._core import Matrix, diag, eye, hessian, jacobian, ones, wronskian, zeros  # noqa: E402

MutableDenseMatrix = Matrix
DenseMatrix = Matrix
class _ImmutableMatrixMeta(type):
    def __instancecheck__(cls, obj):
        return isinstance(obj, Matrix) and obj._immutable

    def __call__(cls, *args):
        return Matrix(*args).as_immutable()


class ImmutableMatrix(metaclass=_ImmutableMatrixMeta):
    """Immutable dense matrix, backed by the native matrix implementation."""

    @staticmethod
    def eye(*args, **kwargs):
        return eye(*args, **kwargs).as_immutable()

    @staticmethod
    def zeros(*args, **kwargs):
        return zeros(*args, **kwargs).as_immutable()

    @staticmethod
    def ones(*args, **kwargs):
        return ones(*args, **kwargs).as_immutable()

    @staticmethod
    def diag(*args, **kwargs):
        return diag(*args, **kwargs).as_immutable()

    @staticmethod
    def hstack(*args):
        return Matrix.hstack(*args).as_immutable()

    @staticmethod
    def vstack(*args):
        return Matrix.vstack(*args).as_immutable()


ImmutableDenseMatrix = ImmutableMatrix


class MatrixError(Exception):
    """Base class for matrix errors (as in SymPy)."""


class ShapeError(ValueError, MatrixError):
    """Wrong matrix shape."""


class NonSquareMatrixError(ShapeError):
    """The matrix must be square."""


class NonInvertibleMatrixError(ValueError, MatrixError):
    """The matrix is singular."""


class NonPositiveDefiniteMatrixError(ValueError, MatrixError):
    """The matrix is not positive-definite."""


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
        return False

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
# --- logic ---
# Boolean algebra (And, Or, Not, ...), SAT, sets and solveset.

from ._core import (  # noqa: E402
    And, Equivalent, ITE, Implies, Nand, Nor, Not, Or, POSform, SOPform, Xnor, Xor,
    ANFform, anf_coeffs, is_anf, to_anf,
    distribute_and_over_or, distribute_or_over_and, eliminate_implications,
    is_cnf, is_dnf, is_literal, is_nnf, simplify_logic, to_cnf, to_dnf, to_nnf,
)
from ._core import satisfiable as _satisfiable, truth_table as _truth_table  # noqa: E402


def satisfiable(expr, algorithm=None, all_models=False, minimal=False, use_lra_theory=False):
    """A model of ``expr`` (dict) or ``False`` (CDCL SAT solver)."""
    return _satisfiable(expr, algorithm, all_models, minimal, use_lra_theory)


def truth_table(expr, variables, input=True):
    """Iterate over the rows of the truth table of ``expr``."""
    return iter(_truth_table(expr, list(variables), input))


def _lg_and(a, b):
    if _core._is_set(a) or _core._is_set(b):
        return Intersection(a, b)
    return And(a, b)


def _lg_or(a, b):
    if _core._is_set(a) or _core._is_set(b):
        return Union(a, b)
    return Or(a, b)


Expr.__and__ = _lg_and
Expr.__rand__ = lambda a, b: _lg_and(b, a)
Expr.__or__ = _lg_or
Expr.__ror__ = lambda a, b: _lg_or(b, a)
Expr.__invert__ = lambda a: Not(a)
Expr.__xor__ = lambda a, b: Xor(a, b)
Expr.__rxor__ = lambda a, b: Xor(b, a)
Expr.__rshift__ = lambda a, b: Implies(a, b)
Expr.__rrshift__ = lambda a, b: Implies(b, a)
Expr.__lshift__ = lambda a, b: Implies(b, a)
Expr.__rlshift__ = lambda a, b: Implies(a, b)


# sets ------------------------------------------------------------------------

_SET_CLASSES = {}


def _wrap_set(e):
    """Wrap a core expression into the matching Set class (if it is a set)."""
    if isinstance(e, Expr):
        cls = _SET_CLASSES.get(e.class_name)
        if cls is not None:
            obj = object.__new__(cls)
            obj._e = e
            return obj
    return e


def _unwrap(o):
    return o._e if isinstance(o, Set) else o


class Set:
    """Base class of sets (wraps a core expression)."""

    __slots__ = ("_e",)
    is_Interval = False
    is_FiniteSet = False
    is_Union = False
    is_Intersection = False
    is_Complement = False
    is_ProductSet = False
    is_EmptySet = False
    is_UniversalSet = False

    def __new__(cls, *args):
        raise TypeError("use a concrete set class")

    # conversions -----------------------------------------------------------
    def _symox_(self):
        return self._e

    def to_sympy(self):
        return self._e.to_sympy()

    def _sympy_(self):
        return self._e.to_sympy()

    def __str__(self):
        return str(self._e)

    __repr__ = __str__

    def srepr(self):
        return self._e.srepr()

    def __eq__(self, other):
        try:
            return self._e == sympify(_unwrap(other))
        except Exception:  # noqa: BLE001
            return False

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self._e)

    def __reduce__(self):
        return (_wrap_set, (self._e,))

    # structure ---------------------------------------------------------------
    @property
    def args(self):
        return tuple(_wrap_set(a) for a in self._e.args)

    @property
    def func(self):
        return type(self)

    @property
    def free_symbols(self):
        return self._e.free_symbols

    @property
    def class_name(self):
        return self._e.class_name

    def subs(self, *args, **kwargs):
        return _wrap_set(self._e.subs(*args, **kwargs))

    def xreplace(self, rule):
        return _wrap_set(self._e.xreplace(rule))

    def evalf(self, *args, **kwargs):
        return _wrap_set(self._e.evalf(*args, **kwargs))

    n = evalf

    def has(self, *patterns):
        return self._e.has(*patterns)

    # set algebra -----------------------------------------------------------
    def union(self, other):
        return Union(self, other)

    def intersect(self, other):
        return Intersection(self, other)

    intersection = intersect

    def complement(self, universe):
        return Complement(universe, self)

    def symmetric_difference(self, other):
        return SymmetricDifference(self, other)

    def __or__(self, other):
        return Union(self, other)

    __ror__ = __or__
    __add__ = __or__
    __radd__ = __or__

    def __and__(self, other):
        return Intersection(self, other)

    __rand__ = __and__

    def __sub__(self, other):
        return Complement(self, other)

    def __rsub__(self, other):
        return Complement(other, self)

    def __xor__(self, other):
        return SymmetricDifference(self, other)

    def __mul__(self, other):
        return ProductSet(self, other)

    def __pow__(self, n):
        n = sympify(n)
        if not n.is_Integer:
            raise ValueError("set powers require a nonnegative integer")
        n = int(n)
        if n < 0:
            raise ValueError("%s: Exponent must be a positive Integer" % n)
        return ProductSet(*([self] * n))

    # membership --------------------------------------------------------------
    def contains(self, other):
        return _core._set_contains(self._e, _unwrap(other))

    def __contains__(self, other):
        r = self.contains(other)
        if r == true:
            return True
        if r == false:
            return False
        raise TypeError("did not evaluate to a bool: %s" % r)

    def as_relational(self, symbol):
        return _core._set_as_relational(self._e, symbol)

    def is_subset(self, other):
        return _core._set_is_subset(self._e, _unwrap(other))

    def issubset(self, other):
        return self.is_subset(other)

    def is_superset(self, other):
        return _core._set_is_subset(_unwrap(other), self._e)

    issuperset = is_superset

    def is_proper_subset(self, other):
        if self == other:
            return False
        return self.is_subset(other)

    def is_proper_superset(self, other):
        if self == other:
            return False
        return self.is_superset(other)

    def isdisjoint(self, other):
        return Intersection(self, other) == EmptySet

    is_disjoint = isdisjoint

    # properties ----------------------------------------------------------------
    @property
    def inf(self):
        return _core._set_inf(self._e)

    @property
    def sup(self):
        return _core._set_sup(self._e)

    @property
    def measure(self):
        return _core._set_measure(self._e)

    @property
    def boundary(self):
        return _wrap_set(_core._set_boundary(self._e))

    @property
    def interior(self):
        return _wrap_set(_core._set_interior(self._e))

    @property
    def closure(self):
        return _wrap_set(_core._set_closure(self._e))

    @property
    def is_empty(self):
        if self.is_EmptySet:
            return True
        if isinstance(self, Range):
            if self.start == self.stop:
                return True
            if all(v.is_Integer for v in self._e.args):
                return False
        if isinstance(self, (Interval, FiniteSet)) or self.class_name in (
                "Reals", "Integers", "Naturals", "Naturals0", "Complexes", "Rationals", "UniversalSet"):
            return False
        return None

    @property
    def is_finite_set(self):
        if isinstance(self, FiniteSet) or self.is_EmptySet:
            return True
        if isinstance(self, Interval) or self.class_name in (
                "Reals", "Integers", "Naturals", "Naturals0", "Complexes", "Rationals"):
            return False
        return None

    @property
    def is_iterable(self):
        return isinstance(self, (FiniteSet, Range)) or self.class_name in ("Integers", "Naturals", "Naturals0")

    @property
    def is_open(self):
        return None

    @property
    def is_closed(self):
        return None

    def __bool__(self):
        return self.is_empty is not True


class Interval(Set):
    """``Interval(start, end, left_open=False, right_open=False)``."""

    __slots__ = ()
    is_Interval = True

    def __new__(cls, start, end, left_open=False, right_open=False):
        return _wrap_set(_core.Interval(_unwrap(start), _unwrap(end), bool(left_open), bool(right_open)))

    @classmethod
    def open(cls, a, b):
        return cls(a, b, True, True)

    @classmethod
    def Lopen(cls, a, b):
        return cls(a, b, True, False)

    @classmethod
    def Ropen(cls, a, b):
        return cls(a, b, False, True)

    @property
    def start(self):
        return self._e.args[0]

    left = start

    @property
    def end(self):
        return self._e.args[1]

    right = end

    @property
    def left_open(self):
        return self._e.args[2]

    @property
    def right_open(self):
        return self._e.args[3]

    @property
    def is_left_unbounded(self):
        return self.start == -oo

    @property
    def is_right_unbounded(self):
        return self.end == oo

    @property
    def is_open(self):
        return bool(self.left_open) and bool(self.right_open)

    @property
    def is_closed(self):
        return (self.start == -oo or not bool(self.left_open)) and (self.end == oo or not bool(self.right_open))


class FiniteSet(Set):
    """``FiniteSet(*elements)``, printed ``{a, b, ...}``."""

    __slots__ = ()
    is_FiniteSet = True

    def __new__(cls, *args):
        return _wrap_set(_core.FiniteSet(*[_unwrap(a) for a in args]))

    def __iter__(self):
        return iter(self.args)

    def __len__(self):
        return len(self._e.args)

    @property
    def is_open(self):
        return False

    @property
    def is_closed(self):
        return True


class Union(Set):
    __slots__ = ()
    is_Union = True

    def __new__(cls, *args):
        if len(args) == 1 and not isinstance(args[0], (Set, Expr)):
            args = tuple(args[0])
        return _wrap_set(_core.Union(*[_unwrap(a) for a in args]))


class Intersection(Set):
    __slots__ = ()
    is_Intersection = True

    def __new__(cls, *args):
        if len(args) == 1 and not isinstance(args[0], (Set, Expr)):
            args = tuple(args[0])
        return _wrap_set(_core.Intersection(*[_unwrap(a) for a in args]))


class Complement(Set):
    __slots__ = ()
    is_Complement = True

    def __new__(cls, a, b):
        return _wrap_set(_core.Complement(_unwrap(a), _unwrap(b)))


def SymmetricDifference(a, b):
    """``(a - b) | (b - a)``."""
    return Union(Complement(a, b), Complement(b, a))


class ProductSet(Set):
    __slots__ = ()
    is_ProductSet = True

    def __new__(cls, *args):
        return _wrap_set(_core.ProductSet(*[_unwrap(a) for a in args]))

    @property
    def sets(self):
        return self.args


class Lambda:
    """``Lambda(var, expr)`` (single variable)."""

    __slots__ = ("_e",)

    def __init__(self, var, expr):
        self._e = _core.Lambda(var, expr)

    def _symox_(self):
        return self._e

    @property
    def free_symbols(self):
        return self._e.free_symbols

    @property
    def variables(self):
        return (self._e.args[0],)

    @property
    def expr(self):
        return self._e.args[1]

    def __call__(self, x):
        return _core._lambda_call(self._e, x)

    def __str__(self):
        return str(self._e)

    __repr__ = __str__

    def __eq__(self, other):
        return isinstance(other, Lambda) and self._e == other._e

    def __hash__(self):
        return hash(self._e)


class ImageSet(Set):
    __slots__ = ()

    def __new__(cls, flambda, *sets):
        base = sets[0] if len(sets) == 1 else ProductSet(*sets)
        return _wrap_set(_core.ImageSet(_unwrap(flambda), _unwrap(base)))

    @property
    def lamda(self):
        lam = object.__new__(Lambda)
        lam._e = self._e.args[0]
        return lam

    @property
    def base_set(self):
        return _wrap_set(self._e.args[1])


class ConditionSet(Set):
    __slots__ = ()

    def __new__(cls, sym, condition, base_set=None):
        b = None if base_set is None else _unwrap(base_set)
        return _wrap_set(_core.ConditionSet(sym, condition, b))

    @property
    def sym(self):
        return self._e.args[0]

    @property
    def condition(self):
        return self._e.args[1]

    @property
    def base_set(self):
        return _wrap_set(self._e.args[2])


class Range(Set):
    __slots__ = ()

    def __new__(cls, *args):
        return _wrap_set(_core.Range(*args))

    @property
    def start(self):
        return self._e.args[0]

    @property
    def stop(self):
        return self._e.args[1]

    @property
    def step(self):
        return self._e.args[2]

    def __iter__(self):
        a, b, s = (int(v) for v in self._e.args)
        return (Integer(i) for i in range(a, b, s))

    def __len__(self):
        a, b, s = (int(v) for v in self._e.args)
        return len(range(a, b, s))


class _NamedSet(Set):
    __slots__ = ()


def _named(name, **flags):
    cls = type(name + "Set", (_NamedSet,), {"__slots__": (), **flags})
    _SET_CLASSES[name] = cls
    return _wrap_set(_core._named_set(name))


for _c in (Interval, FiniteSet, Union, Intersection, Complement, ProductSet, ImageSet,
           ConditionSet, Range):
    _SET_CLASSES[_c.__name__] = _c

EmptySet = _named("EmptySet", is_EmptySet=True,
                  __iter__=lambda self: iter(()), __len__=lambda self: 0)
UniversalSet = _named("UniversalSet", is_UniversalSet=True)
Reals = _named("Reals")
Integers = _named("Integers")
Naturals = _named("Naturals")
Naturals0 = _named("Naturals0")
Complexes = _named("Complexes")
Rationals = _named("Rationals")


def Contains(x, s):
    """Membership ``x in s`` as a boolean expression."""
    return _core.Contains(_unwrap(x), _unwrap(s))


def solveset(f, symbol=None, domain=None):
    """Solve ``f == 0`` (or an (in)equality) for ``symbol`` as a set."""
    return _wrap_set(_core.solveset(f, symbol, None if domain is None else _unwrap(domain)))


def solveset_real(f, symbol):
    return solveset(f, symbol, Reals)


def solveset_complex(f, symbol):
    return solveset(f, symbol, Complexes)


def solve_univariate_inequality(expr, gen, relational=True, domain=None, continuous=False):
    r = _core.solve_univariate_inequality(
        expr, gen, relational, None if domain is None else _unwrap(domain), continuous)
    return _wrap_set(r)


def reduce_inequalities(inequalities, symbols=[]):  # noqa: B006
    if isinstance(inequalities, (list, tuple, set)):
        inequalities = [_unwrap(i) for i in inequalities]
    syms = symbols if symbols else None
    return _core.reduce_inequalities(inequalities, syms)


# --- end logic ---

# --- poly ---
from ._core import (  # noqa: E402
    Poly,
    GroebnerBasis, groebner, is_zero_dimensional,
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
from .printing.pytorch import PyTorchPrinter, TorchPrinter, torchcode  # noqa: E402
from .codegen.pytorch import torchcode_function, torch_lambdify  # noqa: E402
from .printing.jax import JaxPrinter, JAXPrinter, jaxcode  # noqa: E402
from .codegen.jax import jaxcode_function, jax_lambdify  # noqa: E402
from .printing.slang import SlangCodePrinter, SlangPrinter, slangcode  # noqa: E402
from .codegen.slang import slangcode_function, slang_kernel, slang_lambdify  # noqa: E402
from .printing.sycl import SYCLCodePrinter, AdaptiveCppCodePrinter, syclcode, adaptivecppcode  # noqa: E402
from .codegen.sycl import (  # noqa: E402
    syclcode_function, sycl_kernel, adaptivecppcode_function, adaptivecpp_kernel,
)
from .codegen.cuda import CudaError, CudaKernel, compile_cuda, cuda_lambdify  # noqa: E402
from ._sparse import SparseMatrix, MutableSparseMatrix, ImmutableSparseMatrix  # noqa: E402

_native_lambdify = lambdify

_native_N = N


def N(expr, n=15, **kwargs):
    """Evaluate to ``n`` decimal digits using the native precision engine.

    Supports expressions and matrices. Decimal precision is 1 through 10000;
    unsupported evaluation options raise an explicit error.
    """
    if isinstance(expr, _core.Matrix) or getattr(expr, "_symox_sparse", False):
        return expr.evalf(n, **kwargs)
    return _native_N(expr, n, **kwargs)


def lambdify(args, expr, modules=None, **kwargs):
    """Build a numerical callable with native NumPy, PyTorch, JAX or Slang.

    PyTorch supports ``cse``, ``compile``, ``compile_options`` and
    ``user_functions``; JAX supports ``cse``, ``jit``, ``jit_options`` and
    ``user_functions``. Slang uses the optional SlangPy runtime and accepts
    ``device``, ``backend``, ``dtype``, ``cse`` and ``return_type`` settings.
    The default and NumPy paths use the native evaluator.
    """
    candidates = modules if isinstance(modules, (list, tuple)) else [modules]
    torch_modules = [module for module in candidates
                     if (isinstance(module, str) and module in ("torch", "pytorch"))
                     or getattr(module, "__name__", None) == "torch"]
    if torch_modules:
        if len(candidates) != 1:
            raise NotImplementedError("PyTorch lambdify requires a single module; use user_functions for overrides")
        return torch_lambdify(args, expr, **kwargs)
    jax_modules = [module for module in candidates
                   if (isinstance(module, str) and module == "jax")
                   or getattr(module, "__name__", None) == "jax"]
    if jax_modules:
        if len(candidates) != 1:
            raise NotImplementedError("JAX lambdify requires a single module; use user_functions for overrides")
        return jax_lambdify(args, expr, **kwargs)
    slang_modules = [module for module in candidates
                     if (isinstance(module, str) and module in ("slang", "slangpy"))
                     or getattr(module, "__name__", None) == "slangpy"]
    if slang_modules:
        if len(candidates) != 1:
            raise NotImplementedError("Slang lambdify requires a single module; use user_functions for overrides")
        return slang_lambdify(args, expr, **kwargs)
    if any(isinstance(module, str) and module in ("cuda", "nvrtc") for module in candidates):
        if len(candidates) != 1:
            raise NotImplementedError("CUDA lambdify requires a single module")
        return cuda_lambdify(args, expr, **kwargs)
    return _native_lambdify(args, expr, modules=modules, **kwargs)

__all__ = [n for n in dir() if not n.startswith("_")]
