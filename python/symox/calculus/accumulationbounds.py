"""Accumulation bounds: the real interval of the values an expression keeps
taking, as sympy/calculus/accumulationbounds.py (SymPy 1.14).

``AccumBounds(min, max)`` is an expression (a host node of the native tree),
so it appears inside sums, products and powers, prints, compares and
substitutes like the reference class. Its arithmetic, powers, comparisons,
intersection and union are native (``symox::limits::accum``, a line by line
port of the reference methods); the methods here only admit the operands.
"""
from functools import wraps

from .._expr_facade import Expr
from ..multipledispatch import dispatch


def _sympify(value):
    from .. import sympify
    return sympify(value, strict=True)


def _sympifyit(arg, retval=None):
    # sympy.core.decorators._sympifyit (the symox copy imports the public
    # namespace, which is not complete when this module loads).
    def deco(func):
        @wraps(func)
        def wrapper(a, b):
            from ..core.sympify import SympifyError
            try:
                if not hasattr(b, '_op_priority'):
                    b = _sympify(b)
                return func(a, b)
            except SympifyError:
                return retval
        return wrapper
    return deco


def _names():
    from .. import S, FiniteSet
    return S, FiniteSet


# Whether any bounds were built in this process: ``subs`` only looks for
# bounds to re-evaluate once there can be some.
_CREATED = [False]


def _mark_created():
    """Enable the bounds coefficient folding of the Python ``Add`` and
    ``Mul`` constructors once bounds exist (the native computations call it
    when they first return bounds)."""
    if not _CREATED[0]:
        from .. import _compat
        _compat._accum_fold[0] = _fold_coefficient
        _CREATED[0] = True


def _core():
    from .. import _core
    return _core


class AccumulationBounds(Expr):
    r"""An accumulation bounds: the closed real interval ``[min, max]`` of
    the values an expression accumulates at (see the reference docstring)."""

    _symox_user_hooks = True
    __module__ = 'sympy.calculus.accumulationbounds'
    is_extended_real = True
    is_number = False

    def __new__(cls, min, max):
        # The validation and the node are native (symox::limits::accum).
        value = _core()._accum_new(_sympify(min), _sympify(max))
        _mark_created()
        return value

    # setting the operation priority
    _op_priority = 11.0

    def _eval_is_real(self):
        if self.min.is_real and self.max.is_real:
            return True

    def _sympystr(self, printer):
        return "AccumBounds(%s, %s)" % (printer._print(self.min),
                                        printer._print(self.max))

    def to_sympy(self):
        from sympy import AccumBounds
        return AccumBounds(self.min.to_sympy(), self.max.to_sympy())

    @property
    def min(self):
        return self.args[0]

    @property
    def max(self):
        return self.args[1]

    @property
    def delta(self):
        return self.max - self.min

    @property
    def mid(self):
        return (self.min + self.max) / 2

    @_sympifyit('other', NotImplemented)
    def _eval_power(self, other):
        return self.__pow__(other)

    # The arithmetic is native (symox::limits::accum); the methods keep the
    # reference admission of their operands.
    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        if isinstance(other, Expr):
            return _core()._accum_op('add', self, other)
        return NotImplemented

    __radd__ = __add__

    def __neg__(self):
        return _core()._accum_unary('neg', self)

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        if isinstance(other, Expr):
            return _core()._accum_op('sub', self, other)
        return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __rsub__(self, other):
        return _core()._accum_op('rsub', self, other)

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        if self.args == (_names()[0].NegativeInfinity, _names()[0].Infinity):
            return self
        if isinstance(other, Expr):
            return _core()._accum_op('mul', self, other)
        return NotImplemented

    __rmul__ = __mul__

    @_sympifyit('other', NotImplemented)
    def __truediv__(self, other):
        if isinstance(other, Expr):
            return _core()._accum_op('truediv', self, other)
        return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __rtruediv__(self, other):
        if isinstance(other, Expr):
            return _core()._accum_op('rtruediv', self, other)
        return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __pow__(self, other):
        if isinstance(other, Expr):
            return _core()._accum_op('pow', self, other)
        return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __rpow__(self, other):
        return _core()._accum_op('rpow', self, other)

    def __abs__(self):
        return _core()._accum_unary('abs', self)

    def __contains__(self, other):
        """``True`` if ``other`` is contained in self, ``False`` if not,
        otherwise TypeError is raised."""
        other = _sympify(other)
        value = _core()._accum_contains(self, other)
        # The reference returns a Python bool for an infinite element and the
        # evaluated ``And`` (a BooleanAtom) otherwise.
        S = _names()[0]
        if other in (S.Infinity, S.NegativeInfinity):
            return value
        return S.true if value else S.false

    def intersection(self, other):
        """The intersection of ``self`` and ``other`` (an AccumulationBounds
        or a FiniteSet)."""
        S, FiniteSet = _names()
        if not isinstance(other, (AccumBounds, FiniteSet)):
            raise TypeError(
                "Input must be AccumulationBounds or FiniteSet object")
        if isinstance(other, FiniteSet):
            members = _core()._accum_members(self, list(other))
            return FiniteSet(*members) if members else S.EmptySet
        found, value = _core()._accum_intersection(self, other)
        if not found:
            return None
        return S.EmptySet if value is None else value

    @property
    def contains(self):
        # The reference class has no ``contains`` (it is not a Set).
        raise AttributeError("'AccumulationBounds' object has no attribute 'contains'")

    def union(self, other):
        # The reference method (not a true union of disjoint bounds).
        if not isinstance(other, AccumBounds):
            raise TypeError(
                "Input must be AccumulationBounds or FiniteSet object")
        return _core()._accum_union(self, other)


@dispatch(AccumulationBounds, AccumulationBounds)
def _eval_is_ge(lhs, rhs):  # noqa: F811
    return _core()._accum_is_ge(lhs, rhs)


@dispatch(AccumulationBounds, Expr)  # noqa: F811
def _eval_is_ge(lhs, rhs):  # noqa: F811
    return _core()._accum_is_ge(lhs, rhs)


@dispatch(Expr, AccumulationBounds)  # noqa: F811
def _eval_is_ge(lhs, rhs):  # noqa: F811
    return _core()._accum_is_ge(lhs, rhs)


# setting an alias for AccumulationBounds
AccumBounds = AccumulationBounds


def _fold_coefficient(op, args):
    """SymPy ``Add.flatten`` and ``Mul.flatten`` for terms that include
    bounds: the Numbers and the bounds fold into the coefficient through the
    bounds arithmetic (natively), and the coefficient leads the arguments.
    ``None`` leaves the construction to the native engine."""
    nested = 'is_Add' if op == 'Add' else 'is_Mul'
    for arg in args:
        if type(arg) is AccumulationBounds or (getattr(arg, nested, False)
                and any(type(a) is AccumulationBounds for a in arg.args)):
            break
        if op == 'Add' and getattr(arg, 'is_Mul', False) \
                and type(arg.args[0]) is AccumulationBounds:
            break
    else:
        return None
    from .. import sympify
    return _core()._accum_fold(op, [sympify(arg) for arg in args])


_core()._register_accum_bounds(AccumulationBounds._symox_head)
__all__ = ['AccumBounds', 'AccumulationBounds']
