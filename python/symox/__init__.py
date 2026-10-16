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


# --- logic ---
# Boolean algebra (And, Or, Not, ...), SAT, sets and solveset.

from ._core import (  # noqa: E402
    And, Equivalent, ITE, Implies, Nand, Nor, Not, Or, POSform, SOPform, Xnor, Xor,
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
        n = int(n)
        if n < 0:
            raise ValueError("%s: Exponent must be a positive Integer" % n)
        if n == 0:
            raise NotImplementedError("the empty product set is not supported")
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
        return True


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
        return not bool(self.left_open) and not bool(self.right_open)


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

EmptySet = _named("EmptySet", is_EmptySet=True)
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


class _SingletonRegistry:
    """``S``: ``S(obj)`` is ``sympify``; ``S.Reals``, ``S.Half``, ... are the
    singletons."""

    __slots__ = ()

    def __call__(self, *args, **kwargs):
        return sympify(*args, **kwargs)

    def __getattr__(self, name):
        table = {
            "EmptySet": EmptySet, "UniversalSet": UniversalSet, "Reals": Reals,
            "Integers": Integers, "Naturals": Naturals, "Naturals0": Naturals0,
            "Complexes": Complexes, "Rationals": Rationals,
            "true": true, "false": false, "Zero": Integer(0), "One": Integer(1),
            "NegativeOne": Integer(-1), "Half": Rational(1, 2), "Infinity": oo,
            "NegativeInfinity": -oo, "ComplexInfinity": zoo, "NaN": nan, "Pi": pi,
            "Exp1": E, "ImaginaryUnit": I, "EulerGamma": EulerGamma,
        }
        try:
            return table[name]
        except KeyError:
            raise AttributeError(name) from None


S = _SingletonRegistry()
# --- end logic ---


__all__ = [n for n in dir() if not n.startswith("_")]
