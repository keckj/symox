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
    AlgebraicNumber,
    CompiledFunction,
    Derivative,
    Expr,
    Float,
    Integer,
    IsomorphismFailed,
    NotInvertible,
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
    to_number_field,
)

from ._core import Domain, DomainElement, QQ, ZZ, GF, DomainError, CoercionFailed, ExactQuotientFailed, OptionError, FlagError, UnificationFailed

from ._basic import Basic
from . import _node_types
from ._compat import _ExprType
type.__setattr__(Basic, '_sage_', _node_types._sage_placeholder)
S = sympify


class _CallableExpr(Expr):
    """Native lambda and predicate nodes own the callable protocol."""
    __slots__ = ()

    def __call__(self, *args):
        return self._call_symbolic(*args)

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


def _matrix_part(name, M):
    """re, im, conjugate and Abs of a concrete matrix, entrywise as in SymPy."""
    if name in ("re", "im"):
        value = M.as_real_imag()[0 if name == "re" else 1]
    elif name == "conjugate":
        value = M.conjugate()
    else:
        value = M.applyfunc(Abs)
    if type(M).__name__ in ("MutableSparseMatrix", "ImmutableSparseMatrix"):
        return ImmutableSparseMatrix(value)
    return ImmutableMatrix(value)


class _BuiltinFunctionHead:
    """Stable public callable identity and symbolic spelling for native heads."""
    def __init__(self, name, factory):
        import inspect
        self.__name__ = self.__qualname__ = self._symox_head_name = name
        self.__doc__ = factory.__doc__
        self.__module__ = factory.__module__
        self.__signature__ = inspect.signature(factory)
        self._factory = factory

    def __call__(self, *args, **kwargs):
        return self._factory(*args, **kwargs)

    def __getattr__(self, name):
        if name == '_factory':
            raise AttributeError(name)
        # SymPy heads are classes: attributes patched onto them or onto a
        # base class such as Function are visible on the head itself.
        try:
            return _node_types.class_attribute(self.__name__, name)
        except AttributeError:
            return getattr(self._factory, name)

    def __setattr__(self, name, value):
        object.__setattr__(self, name, value)
        if '_factory' in self.__dict__:
            _node_types.forward_setattr(self, self.__name__, name, value)

    def __delattr__(self, name):
        object.__delattr__(self, name)
        _node_types.forward_delattr(self, self.__name__, name)

    def __repr__(self):
        return self.__name__

    __str__ = __repr__

    def __instancecheck__(self, value):
        # SymPy function heads are classes: isinstance(f(x), f).
        return isinstance(value, _core.Expr) and value.func is self

    def __reduce__(self):
        return _public_function_head, (self.__name__,)


def _public_function_head(name):
    return globals()[name]


def _make(name):
    from ._function_signatures import signature as _head_signature,bind as _bind_head
    _contract=_head_signature(name)
    if _contract is not None:
        _fields=list(_contract.parameters.values())
        _required=sum(field.default is field.empty and field.kind!=field.VAR_POSITIONAL for field in _fields)
        _maximum=None if any(field.kind==field.VAR_POSITIONAL for field in _fields) else len(_fields)
        _nullable=tuple(index for index,field in enumerate(_fields) if field.default is None)
    def f(*args, evaluate=None, **kwargs):
        if kwargs:
            if _contract is None:raise TypeError(name+' does not accept these keyword arguments')
            args=_bind_head(_contract,args,kwargs)
        elif _contract is not None:
            if len(args)<_required or _maximum is not None and len(args)>_maximum:
                _contract.bind(*args)  # Raise a Python argument error before entering native code.
            if any(index<len(args) and args[index] is None for index in _nullable):
                args=_bind_head(_contract,args,{})
        from ._public_core import _evaluation
        evaluate = _evaluation.get() if evaluate is None else evaluate
        if name in ('Min','Max'):
            converted=[]
            for arg in args:
                value=sympify(arg)
                numeric_nonreal=(getattr(value,'is_number',False) and getattr(value,'is_real',None) is False
                                 and value not in (oo,-oo))
                if value==nan or numeric_nonreal or getattr(value,'is_extended_real',None) is False or getattr(value,'is_Matrix',False):
                    raise ValueError('%s requires extended-real arguments' % name)
                converted.append(value)
            args=tuple(converted)
        if len(args) == 1 and getattr(args[0], "_symox_matrixexpr", False):
            if name == "conjugate":
                return args[0].conjugate()
            raise NotImplementedError(f"{name} of abstract matrices is not implemented")
        if (name == 'exp' and len(args) == 1 and evaluate
                and (isinstance(args[0], _core.Matrix) or getattr(args[0], '_symox_sparse', False))):
            return args[0].exp().as_immutable()
        if (len(args) == 1 and name in ("re", "im", "conjugate", "Abs")
                and not isinstance(args[0], _core.Expr) and hasattr(args[0], "applyfunc")
                and hasattr(args[0], "as_real_imag")):
            return _matrix_part(name, args[0])
        return _core._apply(name, *args) if evaluate else _core._apply(name, *args, evaluate=False)

    f.__name__ = name
    f.__qualname__ = name
    f.__doc__ = f"symbolic {name}"
    f._symox_head_name = name
    if _contract is not None:
        import inspect as _inspect
        _parameters=list(_contract.parameters.values())
        _parameters.append(_inspect.Parameter('evaluate',_inspect.Parameter.KEYWORD_ONLY,default=None))
        f.__signature__=_contract.replace(parameters=_parameters)
    return _BuiltinFunctionHead(name, f)


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


class _FunctionType(type):
    """Metaclass of Function: ``class F(Function)`` defines a native function class."""

    def __new__(metaclass, name, bases, namespace, **kwargs):
        if "Function" in globals() and any(base is Function for base in bases):
            from ._host import function_class
            return function_class(name, bases, namespace, **kwargs)
        return super().__new__(metaclass, name, bases, namespace, **kwargs)

    def __instancecheck__(cls, value):
        # Applications of Function subclasses are Function instances, as in
        # SymPy; undefined heads keep their nominal class.
        if type.__instancecheck__(cls, value):
            return True
        from ._host import is_host_function_class
        return cls is Function and is_host_function_class(type(value))

    def __subclasscheck__(cls, subclass):
        if type.__subclasscheck__(cls, subclass):
            return True
        from ._host import is_host_function_class
        return cls is Function and is_host_function_class(subclass)

    # SymPy's Function is the base class of every applied head. Patches of
    # it live on its nominal node class, so that undefined heads (instances
    # here, classes in SymPy) and applied nodes resolve them through the MRO.
    def __getattr__(cls, name):
        if cls is Function:
            return _node_types.class_attribute('Function', name)
        raise AttributeError(name)

    def __setattr__(cls, name, value):
        if cls is Function and _node_types._forward_target(cls, 'Function', name) is not None:
            _node_types.forward_setattr(cls, 'Function', name, value)
        else:
            type.__setattr__(cls, name, value)

    def __delattr__(cls, name):
        if cls is Function and name not in cls.__dict__:
            if _node_types._forward_target(cls, 'Function', name) is not None:
                if name not in _node_types.node_class('Function').__dict__:
                    raise AttributeError(name)
                _node_types.forward_delattr(cls, 'Function', name)
                return
        type.__delattr__(cls, name)


def _is_function_head(value):
    """Whether ``value`` is an undefined-function head such as Function('f')."""
    return type.__instancecheck__(Function, value)


class Function(metaclass=_FunctionType):
    """Undefined function: ``f = Function('f'); f(x)``."""

    def __init__(self, name: str, *, nargs=None, **assumptions):
        if not isinstance(name, str):
            raise TypeError("Function expects a string name")
        self.name = name
        self.__name__ = name
        self._assumptions = dict(_core._function_facts(assumptions))
        if nargs is None:
            self._nargs = None
        else:
            import operator
            try:
                counts = (operator.index(nargs),)
            except TypeError:
                counts = tuple(operator.index(n) for n in nargs)
            if not counts or any(n < 0 for n in counts):
                raise ValueError("nargs must contain nonnegative integers")
            self._nargs = tuple(sorted(set(counts)))

    @property
    def nargs(self):
        return Naturals0 if self._nargs is None else FiniteSet(*self._nargs)

    def __call__(self, *args, evaluate=True):
        return _core._undefined_function(self.name, args, self._assumptions, self._nargs)

    def __getattr__(self, name):
        # In SymPy an undefined head is the class of its applications, so
        # class attributes such as SageMath's ``_sage_`` resolve on it.
        if _node_types._patchable(name) and 'name' in self.__dict__:
            from ._function_types import application_type
            cls = application_type(self.name, nargs=self._nargs, **self._assumptions)
            if name in cls.__dict__ or _node_types.has_class_attribute('AppliedUndef', name):
                return getattr(cls, name)
        raise AttributeError("'Function' object has no attribute '%s'" % name)

    def __instancecheck__(self, instance):
        return (isinstance(instance, _core.Expr) and instance.is_UndefinedFunction
                and instance.func == self)

    def __repr__(self):
        return self.name

    def __eq__(self, other):
        return (_is_function_head(other) and other.name == self.name
                and other._assumptions == self._assumptions and other._nargs == self._nargs)

    def __hash__(self):
        return hash(("Function", self.name, tuple(sorted(self._assumptions.items())), self._nargs))

    def __reduce__(self):
        return (_restore_function, (self.name, self._assumptions, self._nargs))

    def to_sympy(self):
        import sympy
        options = dict(self._assumptions)
        if self._nargs is not None:
            options["nargs"] = self._nargs
        return sympy.Function(self.name, **options)

    _sympy_ = to_sympy


def _restore_function(name, assumptions, nargs):
    return Function(name, nargs=nargs, **assumptions)


class _CalculusMeta(type):
    def __getattr__(cls,name):
        try:
            return _node_types.anchored_class_attribute(cls,name)
        except AttributeError:
            pass
        from ._constructor_methods import instance_method
        return instance_method(_core.Expr,name)

    def class_key(cls):return (5,0,cls.__name__)

    def __instancecheck__(cls, value):
        return isinstance(value, Expr) and getattr(value, 'is_'+cls.__name__, False)

    __subclasscheck__ = _ExprType.__subclasscheck__
    __setattr__ = _ExprType.__setattr__
    __delattr__ = _ExprType.__delattr__
    __repr__ = _ExprType.__repr__


class Integral(metaclass=_CalculusMeta):
    """Held integral with ordered symbolic integration limits."""
    def __new__(cls, function, *limits, **kwargs):
        if kwargs:
            raise NotImplementedError("Integral constructor options are not implemented")
        return _core._calculus("Integral", function, limits)


class Sum(metaclass=_CalculusMeta):
    """Held sum with ordered integer limits."""
    def __new__(cls, function, *limits, **kwargs):
        if kwargs:
            raise NotImplementedError("Sum constructor options are not implemented")
        return _core._calculus("Sum", function, limits)


class Product(metaclass=_CalculusMeta):
    """Held product with ordered integer limits."""
    def __new__(cls, function, *limits, **kwargs):
        if kwargs:
            raise NotImplementedError("Product constructor options are not implemented")
        return _core._calculus("Product", function, limits)


from ._concrete import ReorderError  # noqa: E402


class _SubsMeta(_CalculusMeta):
    def __instancecheck__(cls, value):
        return isinstance(value, Expr) and value.is_Subs


class Subs(metaclass=_SubsMeta):
    """Held substitution with explicit bound variables and evaluation points."""
    def __new__(cls, expr, variables, point, **assumptions):
        return _core._subs(expr, variables, point)


# ---------------------------------------------------------------------------
# number theory
# ---------------------------------------------------------------------------


from ._core import (  # noqa: E402
    expand_log, expand_mul, expand_power_base, expand_power_exp, expand_trig,
    logcombine, nsimplify, powsimp, trigsimp,
)


def expand_complex(expr, deep=True):
    """Expand ``expr`` into ``re + I*im`` form, applying only the complex hint."""
    return sympify(expr).expand(deep=deep, complex=True, basic=False, log=False, mul=False,
                                power_exp=False, power_base=False, multinomial=False)


from ._core import checkodesol, dsolve, nsolve, solve  # noqa: E402
from ._core import NonlinearError  # noqa: E402
from ._linear import linear_eq_to_matrix, linsolve, solve_linear_system  # noqa: E402
from ._polysolve import solve_poly_system, nonlinsolve  # noqa: E402
from ._core import product, residue  # noqa: E402
from ._core import (  # noqa: E402
    difference_delta, hypersimilar, hypersimp, rsolve, rsolve_hyper, rsolve_poly, rsolve_ratio,
)

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
    # number theory heads (symox.ntheory holds the integer functions)
    "totient", "reduced_totient", "mobius", "divisor_sigma", "udivisor_sigma",
    "legendre_symbol", "jacobi_symbol", "kronecker_symbol", "primenu",
    "primeomega", "primepi", "carmichael", "tribonacci", "genocchi", "andre",
    "motzkin",
]
for _n in _SPECIAL:
    globals()[_n] = _make(_n)
rf = RisingFactorial  # noqa: F821
ff = FallingFactorial  # noqa: F821

# --- elementary completion ---
# Native heads (symox::special::elementary). Mod/Rem raise ZeroDivisionError
# for a zero divisor like SymPy.
_ELEMENTARY = [
    "Mod", "Rem", "sech", "csch", "asec", "acsc", "asech", "acsch", "acoth",
    "sinc", "frac", "KroneckerDelta", "LeviCivita",
]
for _n in _ELEMENTARY:
    globals()[_n] = _make(_n)
Eijk = LeviCivita  # noqa: F821
Catalan = _core.constant("Catalan")
GoldenRatio = _core.constant("GoldenRatio")
TribonacciConstant = _core.constant("TribonacciConstant")


def _neg_base_odd_root(e):
    # SymPy's real_root Transform: -(-b)**(1/q) for b < 0 and odd q (top-down)
    if (e.is_Pow and e.base.is_negative and e.exp.is_Rational
            and e.exp.p == 1 and e.exp.q % 2):
        return -(-e.base)**e.exp
    if not e.args or not isinstance(e, _core.Expr):
        return e
    args = [_neg_base_odd_root(a) for a in e.args]
    if all(a is b or a == b for a, b in zip(args, e.args)):
        return e
    return e.func(*args)


def real_root(arg, n=None, evaluate=None):
    """Real root of ``arg``: ``real_root(-8, 3) == -2`` (SymPy semantics)."""
    if n is not None:
        arg, n = sympify(arg), sympify(n)  # noqa: F821
        return Piecewise(  # noqa: F821
            (root(arg, n), Or(Eq(n, 1), Eq(n, -1))),  # noqa: F821
            (sign(arg)*root(Abs(arg), n), And(Eq(im(arg), 0), Eq(Mod(n, 2), 1))),  # noqa: F821
            (root(arg, n), True))  # noqa: F821
    return _neg_base_odd_root(sympify(arg))  # noqa: F821
# --- end elementary completion ---


def _tuple(items):
    return _core._apply("Tuple", *items)


def hyper(ap, bq, z):
    """Generalised hypergeometric function ``hyper(ap, bq, z)``."""
    return _core._apply("hyper", _tuple(ap), _tuple(bq), z)


hyper = _BuiltinFunctionHead("hyper", hyper)


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


def Eq(lhs, rhs, **options):
    """Construct a symbolic equality, evaluating it when requested."""
    from ._public_core import _evaluation
    evaluate=options.pop('evaluate',None)
    return _core.Rel("==", lhs, rhs, evaluate=_evaluation.get() if evaluate is None else evaluate)


def _ordered_relation(op,lhs,rhs,options):
    from ._public_core import _evaluation
    evaluate=options.pop('evaluate',None)
    state=_evaluation.get() if evaluate is None else evaluate
    if not state and options:
        name=next(iter(options))
        raise TypeError("unexpected keyword argument %r" % name)
    return _core.Rel(op,lhs,rhs,evaluate=state)


def Ne(lhs, rhs, **options):
    """Construct a symbolic inequality between two expressions."""
    return _ordered_relation('!=',lhs,rhs,options)


def Lt(lhs, rhs, **options):
    """Construct a strict less-than relation."""
    return _ordered_relation('<',lhs,rhs,options)


def Le(lhs, rhs, **options):
    """Construct a less-than-or-equal relation."""
    return _ordered_relation('<=',lhs,rhs,options)


def Gt(lhs, rhs, **options):
    """Construct a strict greater-than relation."""
    return _ordered_relation('>',lhs,rhs,options)


def Ge(lhs, rhs, **options):
    """Construct a greater-than-or-equal relation."""
    return _ordered_relation('>=',lhs,rhs,options)


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
    literals = {}
    for escape in (r'\,', r'\:', r'\ '):
        if escape in names:
            marker = '\0'
            while marker in names or marker in literals or marker.isspace() or marker.isalnum() or marker in ',:\\()':
                marker = chr(ord(marker)+1)
            literals[marker] = escape[1:]
            names = names.replace(escape, marker)
    parts = [p for p in _re.split(r"[\s,]+", names.strip()) if p]
    constructor=Symbol if cls is None else cls
    if not callable(constructor):raise TypeError('Symbol constructor must be callable')
    out = []
    for p in parts:
        for n in _expand_range(p):
            for marker, literal in literals.items():
                n = n.replace(marker, literal)
            out.append(constructor(n,**assumptions))
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
    if isinstance(e, Function):
        return e
    if isinstance(e, type) and any(cls.__name__ == "AppliedUndef" for cls in e.__mro__):
        counts = e.nargs
        count = min(map(int, counts)) if getattr(counts, "is_FiniteSet", False) else 0
        return sympify(e(*(Symbol(f"_function_arg_{i}") for i in range(count)))).func
    return sympify(e)


def to_sympy(e):
    """Convert a symox expression to SymPy."""
    return e.to_sympy() if isinstance(e, Function) else sympify(e).to_sympy()


# --- matrices ---
# Symbolic matrices (Rust implementation in symox._core.Matrix).

from ._core import Matrix, diag, eye, hessian, jacobian, ones, wronskian, zeros  # noqa: E402

MutableDenseMatrix = Matrix
DenseMatrix = Matrix
class _ImmutableMatrixMeta(type):
    def __getattr__(cls,name):
        from ._constructor_methods import instance_method
        return instance_method(Matrix,name)

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


class _CharPolyData:
    """Characteristic-polynomial data for coefficient backends still pending."""

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
        if all(coefficient == 0 for coefficient in self._coeffs):
            return S.NegativeInfinity
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
            return self.all_coeffs() == other.all_coeffs() and str(self.domain) == str(other.domain)
        return False

    def __hash__(self):
        return hash((tuple(self._coeffs), str(self._domain)))

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


from ._purepoly import PurePoly  # noqa: E402


def det(matexpr, **kwargs):
    """Determinant of a matrix."""
    M=matexpr
    if getattr(M, "_symox_matrixexpr", False):
        if kwargs:
            raise NotImplementedError("Custom abstract determinant options are not implemented")
        return M.det()
    return Matrix(M).det(**kwargs)


def trace(expr):
    if getattr(expr, "_symox_matrixexpr", False):
        return expr.trace()
    return Matrix(expr).trace()


def adjoint(arg):
    if not hasattr(arg, "adjoint"):
        arg = sympify(arg)
    return arg.adjoint()


def transpose(arg):
    if getattr(arg, "_symox_matrixexpr", False):
        return arg.transpose()
    if isinstance(arg, (_core.Expr, int, float, complex)):
        # the transpose of a commutative scalar is the scalar itself
        return sympify(arg)
    return Matrix(arg).T


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
    if isinstance(o, Set):
        return o._e
    if isinstance(o, (set, frozenset)):
        return _core.FiniteSet(*(_unwrap(item) for item in o))
    if isinstance(o, tuple):
        return _core._apply('Tuple', *(_unwrap(item) for item in o))
    return o


class Set(metaclass=_node_types.WrapperType):
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
    is_Number = is_Symbol = is_Add = is_Mul = is_Pow = is_Function = False
    is_Integer = is_Rational = is_Float = is_Relational = is_Matrix = False

    @property
    def is_Atom(self):
        return self._e.is_Atom

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
        return isinstance(self, (FiniteSet, Range)) or self.class_name in ("Integers", "Naturals", "Naturals0", "Rationals")

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


BadSignatureError = _core.BadSignatureError
BadArgumentsError = _core.BadArgumentsError
_identity_lambda = None


class _LambdaMeta(type):
    def __call__(cls, *args, **kwargs):
        global _identity_lambda
        if cls is IdentityFunction:
            if kwargs:
                raise TypeError("IdentityFunction takes no constructor options")
            value = _core._apply("IdentityFunction", *args)
        else:
            value = _core.Lambda(*args, **kwargs)
        if value.is_IdentityFunction:
            if _identity_lambda is None:
                _identity_lambda = value
            return _identity_lambda
        return value

    def __instancecheck__(cls, value):
        if type(value) is Lambda or type(value) is IdentityFunction:
            return cls is Lambda or value._e.is_IdentityFunction
        return isinstance(value, Expr) and (value.is_IdentityFunction if cls is IdentityFunction else value.is_Lambda)


def _restore_lambda(value):
    """Upgrade legacy Python Lambda wrappers to the native expression."""
    return IdentityFunction() if value.is_IdentityFunction else value


class Lambda(metaclass=_LambdaMeta):
    """A symbolic function with a tuple signature and exact simultaneous calls.

    Newly constructed values are native Expr objects. The historical slot is
    retained so pickles of the original Python wrapper remain readable.
    """
    __slots__ = ("_e",)

    def __getattr__(self, name):
        if name == "_e":
            raise AttributeError(name)
        return getattr(object.__getattribute__(self, "_e"), name)

    def _symox_(self):
        return self._e

    def _sympy_(self):
        return self._e.to_sympy()

    def __call__(self, *args):
        return self._e(*args)

    def __str__(self):
        return str(self._e)

    __repr__ = __str__

    def __eq__(self, other):
        return self._e == _unwrap(other)

    def __hash__(self):
        return hash(self._e)

    def __reduce__(self):
        return (_restore_lambda, (self._e,))


class IdentityFunction(Lambda):
    __slots__ = ()


class ImageSet(Set):
    __slots__ = ()

    def __new__(cls, flambda, *sets):
        base = sets[0] if len(sets) == 1 else ProductSet(*sets)
        return _wrap_set(_core.ImageSet(_unwrap(flambda), _unwrap(base)))

    @property
    def lamda(self):
        return self._e.args[0]

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
    cls = type(_NamedSet)(name, (_NamedSet,), {"__slots__": (), **flags})
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

from ._set_iteration import install as _install_set_iteration
_install_set_iteration(globals())


def Contains(x, s):
    """Membership ``x in s`` as a boolean expression."""
    return _core.Contains(_unwrap(x), _unwrap(s))


def solveset(f, symbol=None, domain=Complexes):
    """Solve ``f == 0`` (or an (in)equality) for ``symbol`` as a set."""
    return _wrap_set(_core.solveset(f, symbol, None if domain is None else _unwrap(domain)))


def solveset_real(f, symbol):
    return solveset(f, symbol, Reals)


def solveset_complex(f, symbol):
    return solveset(f, symbol, Complexes)


def solve_univariate_inequality(expr, gen, relational=True, domain=Reals, continuous=False):
    def absolute_result():
        condition=sympify(expr)
        if not condition.is_Relational or not condition.has(Abs):
            raise NotImplementedError('No reducible absolute-value condition')
        target=S.Reals if domain is None else _wrap_set(_unwrap(domain))
        if target.is_subset(S.Reals) is not True:
            raise NotImplementedError('Absolute-value inequalities require a real domain')
        from .solvers.inequalities import _abs_domain
        return _abs_domain(condition.lhs-condition.rhs,condition.rel_op,gen).intersect(target)

    try:
        r = _core.solve_univariate_inequality(
            expr, gen, False, None if domain is None else _unwrap(domain), continuous)
    except NotImplementedError as native_error:
        try:result=absolute_result()
        except (ValueError,TypeError,NotImplementedError):raise native_error
    else:
        result=_wrap_set(r)
        if isinstance(result,ConditionSet):
            try:result=absolute_result()
            except (ValueError,TypeError,NotImplementedError):pass
    return result.as_relational(gen) if relational else result


def reduce_inequalities(inequalities, symbols=[]):  # noqa: B006
    if isinstance(inequalities, (list, tuple, set)) or hasattr(inequalities,'__next__'):
        inequalities = [_unwrap(i) for i in inequalities]
    if hasattr(symbols,'__next__'):symbols=list(symbols)
    syms = symbols if symbols else None
    conditions=inequalities if isinstance(inequalities,(list,tuple,set)) else [inequalities]
    if any(sympify(condition).has(Abs) for condition in conditions):
        from .solvers.inequalities import _reduce_abs_system
        try:return _reduce_abs_system(inequalities,syms)
        except (ValueError,TypeError,NotImplementedError):pass
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
    minimal_polynomial,
    minpoly,
    primitive_element,
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
    subresultants,
    sturm, count_roots, intervals, refine_root, ground_roots, gff_list, RefinementFailed,
    roots,
    sqf,
    sqf_list,
    sqf_part,
    terms_gcd,
)

together = _core._poly_together


from .polys.polyerrors import (  # noqa: E402
    BasePolynomialError, ComputationFailed, EvaluationFailed, ExtraneousFactors,
    GeneratorsError, GeneratorsNeeded, HeuristicGCDFailed, HomomorphismFailed,
    ModularGCDFailed, MultivariatePolynomialError, NotAlgebraic, NotReversible,
    OperationNotSupported, PolificationFailed, PolynomialDivisionFailed,
    PolynomialError, UnivariatePolynomialError, UnsolvableFactorError,
)


def gcd_list(seq, *gens, **args):
    return gcd(list(seq), *gens, **args)


def lcm_list(seq, *gens, **args):
    return lcm(list(seq), *gens, **args)


def poly(expr, *gens, **args):
    return Poly(expr, *gens, **args)


def primitive(f, *gens, **args):
    """``(content, primitive part)`` of a polynomial."""
    c, p = Poly(f, *gens, **args).primitive()
    return c, p.as_expr()


def _install_poly_methods():
    def as_poly(self,*gens,**kw):
        try:
            return Poly(self,*gens,**kw)
        except PolynomialError:
            return None
        except ValueError as error:
            # The native constructor still reports this conversion failure as
            # ValueError. Do not swallow invalid options or duplicate gens.
            if 'contains an element of the set of generators' in str(error):
                return None
            raise
    methods = {
        "factor": lambda self, *gens, **kw: factor(self, *gens, **kw),
        "cancel": lambda self, *gens, **kw: cancel(self, *gens, **kw),
        "apart": lambda self, x=None, **kw: apart(self, x, **kw),
        "collect": lambda self, syms, **kw: collect(self, syms, **kw),
        "as_poly": as_poly,
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

from .polys import _poly_methods  # noqa: E402
_poly_methods.install()
from .polys import _evaluation as _poly_evaluation  # noqa: E402
_poly_evaluation.install()
from ._purepoly import install as _install_purepoly_methods  # noqa: E402
_install_purepoly_methods()
from .polys.numberfields.galoisgroups import galois_group  # noqa: E402
from .polys.constructor import construct_domain  # noqa: E402
from .polys.polytools import (  # noqa: E402
    LC, LM, LT, all_roots, content, decompose, degree_list, gff, monic,
    nth_power_roots_poly, parallel_poly_from_expr, poly_from_expr, real_roots,
    sqf_norm, total_degree, trunc,
)
from .polys.polyfuncs import (  # noqa: E402
    horner, interpolate, rational_interpolate, symmetrize, viete,
)
from .polys.specialpolys import (  # noqa: E402
    cyclotomic_poly, interpolating_poly, random_poly, swinnerton_dyer_poly,
    symmetric_poly,
)
from .polys.orthopolys import (  # noqa: E402
    chebyshevt_poly, chebyshevu_poly, hermite_poly,
    hermite_prob_poly, jacobi_poly, laguerre_poly, legendre_poly,
)
from .polys.partfrac import apart_list, assemble_partfrac_list  # noqa: E402
from .polys.monomials import Monomial, itermonomials  # noqa: E402
from .polys.orderings import grevlex, grlex, igrevlex, igrlex, ilex, lex  # noqa: E402
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


class Order(metaclass=_ExprType):
    """Landau order term ``O(expr, (x, x0), ...)``."""
    _matches = staticmethod(lambda expr: expr.class_name == "Order")

    def __new__(cls, expr, *args, **kwargs):
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
        if self.is_Lambda:
            if i == 0:
                return self._lambda_body
            if i == 1:
                return self._lambda_variables
        if (self.is_Integral or self.is_Sum or self.is_Product) and i == 1:
            return [limit[0] for limit in self.limits]
        if self.is_Derivative:
            if i == 0:
                return self.args[0]
            if i == 1:
                return tuple(v for v, n in self.variable_count for _ in range(n))
        if self.is_Limit:
            if i == 0:
                return self.args[0]
            if i == 2:
                return self.args[2]
        if self.is_Subs:
            return self.args[0] if i == 0 else tuple(self.args[i].args)
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
Expr.var = property(lambda self: self.args[1] if self.is_Limit else (_ for _ in ()).throw(AttributeError("only limits have var")))
Expr.dir = property(lambda self: str(self.args[3]) if self.is_Limit else (_ for _ in ()).throw(AttributeError("only limits have dir")))

class _LimitMeta(type):
    def __instancecheck__(cls, value):
        return isinstance(value, Expr) and value.is_Limit


class Limit(metaclass=_LimitMeta):
    """Held native limit; ``doit()`` invokes the exact limit engine."""
    def __new__(cls, e, z, z0, dir="+"):
        return _core._held_limit(e, z, z0, str(dir))


def limit(e, z, z0, dir="+"):
    """Limit of ``e`` as ``z -> z0`` (``dir`` is ``'+'``, ``'-'`` or ``'+-'``).

    Returns an unevaluated :class:`Limit` when the limit cannot be computed.
    """
    r = _core._limit(e, z, z0, str(dir))
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
    from ._cse_matrices import is_concrete_matrix,flatten
    single = isinstance(exprs, (Expr, int, float, str)) or is_concrete_matrix(exprs)
    if isinstance(exprs, dict):
        keys = _builtins.list(exprs)
        reps, red = cse([exprs[k] for k in keys], symbols, optimizations, postprocess, order, ignore)
        return reps, dict(zip(keys, red))
    seq = [exprs] if single else exprs
    seq,restore=flatten(seq)
    opt = optimizations not in (False, "none")
    reps, red = _core._cse(seq, opt, ignore or None,
                         order='none' if order=='none' else 'canonical')
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
    red=restore(red)
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
    Add, Mul, Pow, Symbol, Dummy, Wild, Number, NumberSymbol, Integer, Rational, Float, Tuple, UnevaluatedExpr,
    S, count_ops, default_sort_key,
)
from ._compat import (  # noqa: E402
    Atom, AtomicExpr, Relational, Rel, Equality, Unequality, GreaterThan, LessThan,
    StrictGreaterThan, StrictLessThan,
)
from ._compat import _install_relation_classes
_install_relation_classes(globals())
Id = S.IdentityFunction
from .codegen.ast import Indexed, IndexedBase  # noqa: E402
from .printing.str import StrPrinter  # noqa: E402
from .printing.pretty import PrettyPrinter  # noqa: E402
from . import extra  # noqa: E402
from ._sparse import SparseMatrix, MutableSparseMatrix, ImmutableSparseMatrix  # noqa: E402
from . import _matrix_aux, _sparse as _sparse_module  # noqa: E402
from ._matrix_aux import (  # noqa: E402
    GramSchmidt, casoratian, matrix_multiply_elementwise, list2numpy, matrix2numpy, symarray, randMatrix,
    rot_givens, rot_axis1, rot_axis2, rot_axis3, rot_ccw_axis1, rot_ccw_axis2, rot_ccw_axis3, banded,
    MatrixBase, jordan_cell_function as jordan_cell,
)
_matrix_aux.install(ImmutableMatrix, _sparse_module._SparseBase, MutableSparseMatrix)
from . import _matrix_decompositions  # noqa: E402
_matrix_decompositions.install(Matrix, _sparse_module._SparseBase)
from . import _matrix_charpoly  # noqa: E402
_matrix_charpoly.install(Matrix)
MutableMatrix = Matrix

_native_lambdify = lambdify

_native_N = N


def N(expr, n=15, **kwargs):
    """Evaluate to ``n`` decimal digits using the native precision engine.

    Supports expressions and matrices, with precision retained in native numbers.
    Unsupported evaluation options raise an explicit error.
    """
    if getattr(expr, "_symox_matrixexpr", False):
        raise NotImplementedError("Materialize abstract matrices with as_explicit() before numerical evaluation")
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
        return extra.torch_lambdify(args, expr, **kwargs)
    jax_modules = [module for module in candidates
                   if (isinstance(module, str) and module == "jax")
                   or getattr(module, "__name__", None) == "jax"]
    if jax_modules:
        if len(candidates) != 1:
            raise NotImplementedError("JAX lambdify requires a single module; use user_functions for overrides")
        return extra.jax_lambdify(args, expr, **kwargs)
    slang_modules = [module for module in candidates
                     if (isinstance(module, str) and module in ("slang", "slangpy"))
                     or getattr(module, "__name__", None) == "slangpy"]
    if slang_modules:
        if len(candidates) != 1:
            raise NotImplementedError("Slang lambdify requires a single module; use user_functions for overrides")
        return extra.slang_lambdify(args, expr, **kwargs)
    if any(isinstance(module, str) and module in ("cuda", "nvrtc") for module in candidates):
        if len(candidates) != 1:
            raise NotImplementedError("CUDA lambdify requires a single module")
        return extra.cuda_lambdify(args, expr, **kwargs)
    from .utilities.lambdify import _compile_implemented,_compile_structured
    structured = _compile_structured(args,expr,modules,kwargs)
    if structured is not None:
        return structured
    implemented = _compile_implemented(args, expr, modules, kwargs)
    if implemented is not None:
        return implemented
    return _native_lambdify(args, expr, modules=modules, **kwargs)

from ._matrixexpr import (  # noqa: E402
    MatrixExpr, MatrixSymbol, FunctionMatrix, ElementwiseApplyFunction, MatAdd, MatMul, MatPow, Transpose, Adjoint,
    Inverse, Identity, ZeroMatrix, MatrixElement, MatrixConjugate, Trace, Determinant,
    OneMatrix, HadamardProduct, HadamardPower, KroneckerProduct,
    BlockMatrix, BlockDiagMatrix, MatrixSlice, block_collapse, DotProduct,
    DiagMatrix, DiagonalOf, DiagonalMatrix,
    hadamard_product, hadamard_power, kronecker_product,
)
from ._matrix_permutations import PermutationMatrix, MatrixPermute  # noqa: E402
from ._simptools import (  # noqa: E402
    fraction, numer, denom, collect_const, rcollect, radsimp, signsimp, posify, separatevars,
    powsimp, powdenest, factor_terms, gcd_terms, factor_nc, sqrtdenest, combsimp, gammasimp,
    ratsimp, ratsimpmodprime, fu, exptrigsimp, trigsimp, FU, kroneckersimp,
)

# Assumption-dependent queries and rewrites use the native Rust engine.
from .assumptions import (Q, Predicate, AppliedPredicate, AssumptionsContext,
                          global_assumptions, assuming, ask, refine)
Expr.refine = lambda self, assumptions=True: refine(self, assumptions)


def _rewrite_name(r):
    # SymPy dispatches on `_eval_rewrite_as_<name>`; `abs` means `Abs`.
    if r is abs:
        return "Abs"
    if isinstance(r, str):
        return r
    n = getattr(r, "__name__", None)
    return n if n is not None else type(r).__name__


from . import _host  # noqa: E402


def _expr_rewrite(self, *args, deep=True, **hints):
    """Rewrite in terms of ``args[-1]``, restricted to instances of
    ``args[:-1]`` when given (SymPy's ``Basic.rewrite``)."""
    if not args:
        return self
    pattern, rule = args[:-1], args[-1]
    if pattern and isinstance(pattern[0], (list, tuple, set, frozenset)):
        pattern = tuple(pattern[0])
    if self.is_Piecewise and _rewrite_name(rule) == 'ITE' and (
            not pattern or any(_rewrite_name(p) == 'Piecewise' for p in pattern)):
        return _core._boolean_piecewise_to_ite(self)
    names = [_rewrite_name(p) for p in pattern]
    rule_name = _rewrite_name(rule)
    if _host.host_classes_defined() and _core._host_contains(self)[0]:
        if pattern:
            pattern = tuple(p for p in pattern if self.has(p))
            if not pattern:
                return self
        return _host.rewrite(self, pattern, rule, '_eval_rewrite_as_' + rule_name,
                             lambda e, d: _core._rewrite(e, names, rule_name, d), deep, hints)
    return _core._rewrite(self, names, rule_name, deep)


Expr.rewrite = _expr_rewrite

from . import _expr_methods  # noqa: E402
exp_polar = _make('exp_polar')
from .series import (  # noqa: E402
    series, fps, fourier_series, FormalPowerSeries, FourierSeries, FiniteFourierSeries,
    SeqFormula, SeqPer, SeqAdd, SeqMul, EmptySequence, sequence,
    SeqExpr, SeqExprOp, RecursiveSeq,
    approximants, limit_seq, gruntz,
)
_expr_methods.install()

from . import ntheory, discrete  # noqa: E402
from .ntheory import (  # noqa: E402
    Sieve, abundance, binomial_coefficients, binomial_coefficients_list,
    composite, compositepi, continued_fraction, continued_fraction_convergents,
    continued_fraction_iterator, continued_fraction_periodic,
    continued_fraction_reduce, crt, cycle_length, discrete_log, divisor_count,
    divisors, egyptian_fraction, factor_cache, factorint, factorrat, igcd, ilcm,
    integer_log, integer_nthroot, is_abundant, is_amicable, is_carmichael,
    is_deficient, is_mersenne_prime, is_nthpow_residue, is_perfect,
    is_primitive_root, is_quad_residue, isprime, mersenne_prime_exponent,
    mod_inverse, multinomial_coefficients, multiplicity, n_order, nextprime,
    npartitions, nthroot_mod, num_digits, perfect_power, pollard_pm1,
    pollard_rho, prevprime, prime, primefactors, primerange, primitive_root,
    primorial, proper_divisor_count, proper_divisors, quadratic_congruence,
    quadratic_residues, randprime, sieve, sqrt_mod, sqrt_mod_iter, trailing)
from .discrete import (  # noqa: E402
    convolution, covering_product, fft, fwht, ifft, ifwht, intersecting_product,
    intt, inverse_mobius_transform, mobius_transform, ntt)

if hasattr(_core, 'EX'):
    from ._core import EX, EXRAW, RR, CC, ZZ_I, QQ_I
    from .polys.domains._types import (
        IntegerRing, RationalField, FiniteField, RealField, ComplexField,
        GaussianIntegerRing, GaussianRationalField, ExpressionDomain, ExpressionRawDomain,
        PolynomialRing, FractionField, AlgebraicField,
    )
    FF = GF

from .polys.domains import _methods as _domain_methods  # noqa: E402
_domain_methods.install()

from . import _matrixexpr_methods  # noqa: E402
_matrixexpr_methods.install()
from . import _matrixexpr_algebra  # noqa: E402
_matrixexpr_algebra.install()
from . import _algebraic_number_methods  # noqa: E402
_algebraic_number_methods.install()

# Integral transforms (Laplace, Fourier, sine/cosine, Mellin, Hankel).
from ._transforms import (  # noqa: E402
    IntegralTransform, IntegralTransformError, LaplaceTransform, InverseLaplaceTransform,
    FourierTransform, InverseFourierTransform, SineTransform, InverseSineTransform,
    CosineTransform, InverseCosineTransform, HankelTransform, InverseHankelTransform,
    MellinTransform, InverseMellinTransform, laplace_transform, inverse_laplace_transform,
    fourier_transform, inverse_fourier_transform, sine_transform, inverse_sine_transform,
    cosine_transform, inverse_cosine_transform, hankel_transform, inverse_hankel_transform,
    mellin_transform, inverse_mellin_transform, laplace_correspondence, laplace_initial_conds,
)

from ._differentiation import diff, install as _install_differentiation  # noqa: E402
_install_differentiation()
from ._numeric_recognition import nsimplify  # noqa: E402
from ._numeric_recognition import _limit_denominator
_core.Expr.limit_denominator = _limit_denominator
from .solvers.inequalities import (  # noqa: E402
    solve_poly_inequality, solve_rational_inequalities,
    reduce_abs_inequality, reduce_abs_inequalities,
)

from .tensor.array import (  # noqa: E402
    Array, NDimArray, DenseNDimArray, SparseNDimArray, MutableNDimArray,
    ImmutableNDimArray, MutableDenseNDimArray, ImmutableDenseNDimArray,
    MutableSparseNDimArray, ImmutableSparseNDimArray, tensorproduct,
    tensorcontraction, tensordiagonal, permutedims, derive_by_array,
    ArrayComprehension, ArrayComprehensionMap, ArrayKind,
)
from ._kind import Kind,NumberKind,BooleanKind,UndefinedKind
from .calculus.accumulationbounds import AccumBounds,AccumulationBounds
from .tensor.functions import shape  # noqa: E402
from .tensor.array import _adapters as _array_adapters  # noqa: E402
_array_adapters.install(globals())

from .calculus.finite_diff import (  # noqa: E402
    finite_diff_weights, apply_finite_diff, differentiate_finite,
)
from .calculus import finite_diff as _finite_diff  # noqa: E402
_finite_diff.install()
from .calculus import (singularities, stationary_points, periodicity,  # noqa: E402
    is_increasing, is_strictly_increasing, is_decreasing,
    is_strictly_decreasing, is_monotonic, is_convex, minimum, maximum, not_empty_in)

from .utilities import (  # noqa: E402
    ordered, sift, flatten, unflatten, reshape, group, subsets, variations,
    cartes, prefixes, postfixes, rotations, has_dups, has_variety, take, prod,
    dict_merge, topological_sort, preorder_traversal, postorder_traversal,
    bottom_up, use, cacheit, memoize_property,
)
from ._polar import (polar_lift, periodic_argument, principal_branch,
                     unbranched_argument, polarify, unpolarify)
from .utilities.iterables import (numbered_symbols, filter_symbols,
    multiset, multiset_combinations, multiset_permutations,
    multiset_partitions, partitions, ordered_partitions,
    common_prefix, common_suffix, rotate_left, rotate_right,
    interactive_traversal, iterable, is_sequence)

from .plotting import (plot, plot_parametric, plot3d, plot3d_parametric_line,
    plot3d_parametric_surface, plot_contour, plot_implicit, textplot, PlotGrid)

from ._solver_checks import (  # noqa: E402
    checksol, solve_linear, solve_linear_system_LU, solve_undetermined_coeffs,
)

from .solvers.decompogen import decompogen  # noqa: E402
from .solvers.deutils import ode_order, homogeneous_order  # noqa: E402
from ._numeric_float import nfloat  # noqa: E402
from .solvers.ode import classify_ode, dsolve, checkodesol  # noqa: E402
from .solvers.solvers import det_quick,inv_quick
from .solvers.pde import (  # noqa: E402
    pdsolve, classify_pde, checkpdesol, pde_separate, pde_separate_add, pde_separate_mul,
)
from .geometry import (  # noqa: E402
    Point, Point2D, Point3D, Line, Ray, Segment, Line2D, Ray2D, Segment2D,
    Line3D, Ray3D, Segment3D, Plane, Ellipse, Circle, Polygon, Triangle,
    RegularPolygon, Curve, Parabola, rad, deg, are_similar, centroid, convex_hull,
    idiff, intersection, closest_points, farthest_points, GeometryError,
)

from .tensor import _adapters as _indexed_tensor_adapters
_indexed_tensor_adapters.install(globals())


from ._epath import EPath, epath  # noqa: E402
from .polys.rings import ring, xring, vring, sring  # noqa: E402
from .polys.fields import field, xfield, vfield, sfield  # noqa: E402
from ._set_extensions import ComplexRegion, DisjointUnion, PowerSet, imageset  # noqa: E402
from .sets.ordinals import Ordinal, OmegaPower, ord0  # noqa: E402

from .polys.numberfields import field_isomorphism, isolate, round_two, prime_decomp, prime_valuation  # noqa: E402

from ._special_extensions import (
    intlike, multigamma, TupleArg, TupleParametersBase, meijerg, appellf1,
    BesselBase, SphericalBesselBase, SphericalHankelBase, AiryBase, assume_integer_order,
    hn1, hn2, jn_zeros, marcumq, MathieuBase, mathieuc, mathieus, mathieucprime, mathieusprime,
    HyperRep, HyperRep_power1, HyperRep_power2, HyperRep_log1, HyperRep_log2,
    HyperRep_atanh, HyperRep_asin1, HyperRep_asin2, HyperRep_sqrts1, HyperRep_sqrts2,
    HyperRep_cosasin, HyperRep_sinasin,
)
from . import _special_extensions
_special_extensions.install()
from . import functions  # noqa: E402

from .parsing.sympy_parser import parse_expr  # noqa: E402
from .printing import (  # noqa: E402
    TableForm, dotprint, glsl_code, jscode, julia_code, latex, maple_code,
    mathematica_code, mathml, multiline_latex, octave_code, pager_print,
    pprint, pretty, pretty_print, pprint_try_use_unicode, pprint_use_unicode,
    preview, print_fcode, print_glsl, print_gtk, print_jscode, print_latex,
    print_maple_code, print_mathml, print_python, print_rcode, print_tree,
    python, rcode, smtlib_code, sstr, sstrrepr,
)

# Complete the public surface with native-backed implementations.
from . import algebras, concrete, external, interactive, multipledispatch, release, strategies
from .algebras import Quaternion
from .polys.domains._compat import (FF_gmpy,FF_python,GMPYFiniteField,GMPYIntegerRing,
    GMPYRationalField,PythonFiniteField,PythonIntegerRing,PythonRational,
    QQ_gmpy,QQ_python,ZZ_gmpy,ZZ_python)
from ._public_tools import (arity,comp,capture,filldedent,public,threaded,xthreaded,timed,
    vectorize,failing_assumptions,check_assumptions,register_handler,remove_handler,
    bool_map,expand_multinomial)
from ._public_core import (Dict,FunctionClass,WildFunction,PrecisionExhausted,
    SympifyError,evalf,evaluate,seterr)
RealNumber = Float
from ._public_matrix import (DeferredVector,Idx,Permanent,per,blockcut,
    diagonalize_vector,matrix_symbols,get_indices,get_contraction_structure)
from ._public_calculus import (piecewise_fold,piecewise_exclusive,bspline_basis,
    bspline_basis_set,interpolating_spline,line_integrate,euler_equations)
from ._public_special import (betainc,betainc_regularized,erf2inv,Ynm_c,Znm,
    chebyshevt_root,chebyshevu_root,hermite_prob,jacobi_normalized,riemann_xi,
    stieltjes,SingularityFunction,singularityintegrate)
from ._public_roots import RootOf,CRootOf,ComplexRootOf,rootof,RootSum
from ._public_solvers import diophantine,factor_system,solve_triangulated,substitution
from .polys.polyoptions import Options
from .interactive import init_printing,init_session
from .plotting import plot_backends
from ._public_testing import test,doctest
from . import _public_core,_public_matrix,_public_special,_public_roots
_public_core.install(globals())
_public_matrix.install(globals())
_public_roots.install()
_public_special.install(globals())
from . import _combinatorial_rewrites  # noqa: E402
_combinatorial_rewrites.install()
from ._tuple_protocol import install as _install_tuple_protocol
_install_tuple_protocol(Expr)
from ._boolean_operators import install as _install_boolean_operators
_install_boolean_operators()
_native_piecewise = Piecewise

from ._compat import _ExprType as _PiecewiseBaseMeta
class _PiecewiseConstructorMeta(_PiecewiseBaseMeta):
    def __repr__(cls):return cls.__name__
    __str__=__repr__
    def class_key(cls):return (4,0,cls.__name__)

class Piecewise(metaclass=_PiecewiseConstructorMeta):
    """Construct ordered conditional branches with an explicit evaluation mode."""
    _matches=staticmethod(lambda expression:expression.is_Piecewise)
    is_Piecewise=True
    def __new__(cls,*pieces,evaluate=None):
        if not pieces:raise TypeError('At least one (expr, cond) pair expected.')
        state=_public_core._evaluation.get() if evaluate is None else evaluate
        if state:
            return _native_piecewise(*pieces,evaluate=True)
        # Even held Piecewise constructors discard false branches and stop
        # at the first true branch, as required by the branch protocol.
        held=_native_piecewise(*pieces,evaluate=False)
        branches=[]
        for pair in held.args:
            if pair.cond==false:
                continue
            branches.append(pair)
            if pair.cond==true:
                break
        if len(branches)==1 and branches[0].cond==true:
            return branches[0].expr
        return _native_piecewise(*branches,evaluate=False)

import sys as _piecewise_sys
for _piecewise_module_name,_piecewise_module in tuple(_piecewise_sys.modules.items()):
    if (_piecewise_module_name.startswith('symox.') and _piecewise_module is not None
            and _piecewise_module.__dict__.get('Piecewise') is _native_piecewise):
        _piecewise_module.Piecewise=Piecewise


for _context_name in ('sqrt','cbrt','log'):
    def _context_function(name,original):
        def apply(args,evaluate):
            state = _public_core._evaluation.get() if evaluate is None else evaluate
            if state:return original(*args)
            if name in ('sqrt','cbrt'):
                return _core._node('Pow',(args[0],Rational(1,2 if name=='sqrt' else 3)),False)
            if name == 'log' and len(args) == 2:
                return _core._undefined_function('log',args,{},(1,2))
            return _core._apply(name,*args,evaluate=False)
        if name == 'log':
            def f(arg,base=None,*,evaluate=None):
                return apply((arg,) if base is None else (arg,base),evaluate)
            f._symox_head_name = name
            f._nargs = (1,2)
        else:
            def f(arg,evaluate=None):
                return apply((arg,),evaluate)
        f.__name__ = name
        f.__doc__ = original.__doc__ or 'Construct symbolic '+name+'.'
        return _BuiltinFunctionHead(name, f)
    globals()[_context_name] = _context_function(_context_name,globals()[_context_name])
ln = log

def root(arg,n,k=0,evaluate=None):
    n = sympify(n)
    principal = Pow(arg,S.One/n,evaluate=evaluate)
    if k:
        return Mul(principal,S.NegativeOne**(2*k/n),evaluate=evaluate)
    return principal

from ._general_simplify import install as _install_general_simplify
_install_general_simplify(globals())

from ._solve_relations import install as _install_solve_relations
_install_solve_relations(globals())

from .printing.repr import srepr

__all__ = [n for n in dir() if not n.startswith("_")]

# Preload the compatible package before restoring the public callable. A later
# import of symox.simplify must not overwrite the scalar simplification API.
from importlib import import_module as _import_api_module
_simplify_callable = simplify
_import_api_module('.simplify', __name__)
simplify = _simplify_callable

# SymPy exposes its already initialized core namespace from a root import.
# These are the native modules, so downstream type checks retain identity.
_import_api_module('.core', __name__)
_import_api_module('.core.expr', __name__)

from .utilities.cache import CACHE as _public_caches
from ._compat import _cached_symbol
_public_caches.append(_cached_symbol)

from ._function_types import application_type as _function_application_type

from ._distribution_methods import install as _install_distribution_methods
_install_distribution_methods(globals())

from ._matrix_subclasses import install as _install_matrix_subclasses
_install_matrix_subclasses()
