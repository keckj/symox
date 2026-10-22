"""Accumulation bounds: the real interval of the values an expression keeps
taking, as sympy/calculus/accumulationbounds.py (SymPy 1.14).

``AccumBounds(min, max)`` is an expression (a host node of the native tree),
so it appears inside sums, products and powers, prints, compares and
substitutes like the reference class. Its arithmetic, powers, comparisons,
intersection and union follow the reference methods line by line.
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
    from .. import Add, Mul, Pow, S, Min, Max, And, FiniteSet, oo, zoo
    return Add, Mul, Pow, S, Min, Max, And, FiniteSet, oo, zoo


def _extended_real(e):
    """``e.is_extended_real``, completed for the products, sums and
    powers whose realness the kernel leaves unknown but the reference
    derives: ``oo*(r - 1)`` and ``nn**nn`` with a real r and a nonnegative
    nn are extended real."""
    value = e.is_extended_real
    if value is not None or not getattr(e, 'args', ()):
        return value
    if e.is_Mul or e.is_Add:
        return True if all(_extended_real(a) for a in e.args) else None
    if e.is_Pow:
        base, exp = e.args
        if base.is_extended_nonnegative and exp.is_extended_real and base.is_finite \
                and exp.is_finite:
            return True
    return None


# Whether any bounds were built in this process: ``subs`` only looks for
# bounds to re-evaluate once there can be some.
_CREATED = [False]


def _flat_mul(*factors):
    """The reference ``Mul(*factors)``: ``Mul.flatten`` multiplies numbers
    and bounds into one leading coefficient with the bounds arithmetic, which
    the native product does not do for host nodes."""
    Add, Mul, Pow, S, Min, Max, And, FiniteSet, oo, zoo = _names()
    coeff, rest = S.One, []
    for factor in factors:
        for term in Mul.make_args(factor):
            if isinstance(term, AccumulationBounds):
                coeff = term.__mul__(coeff)
            elif term.is_Number:
                coeff = coeff*term
            else:
                rest.append(term)
    if not rest:
        return coeff
    product = Mul(*rest)
    if not isinstance(coeff, AccumulationBounds):
        return coeff*product
    return Mul(coeff, *Mul.make_args(product), evaluate=False)


class AccumulationBounds(Expr):
    r"""An accumulation bounds: the closed real interval ``[min, max]`` of
    the values an expression accumulates at (see the reference docstring)."""

    _symox_user_hooks = True
    __module__ = 'sympy.calculus.accumulationbounds'
    is_extended_real = True
    is_number = False

    def __new__(cls, min, max):
        Add, Mul, Pow, S, Min, Max, And, FiniteSet, oo, zoo = _names()
        min = _sympify(min)
        max = _sympify(max)

        # Only allow real intervals (use symbols with 'is_extended_real=True').
        if not _extended_real(min) or not _extended_real(max):
            raise ValueError("Only real AccumulationBounds are supported")

        if max == min:
            return max

        # Make sure that the created AccumBounds object will be valid.
        if max.is_number and min.is_number:
            bad = max.is_comparable and min.is_comparable and max < min
        else:
            bad = (max - min).is_extended_negative
        if bad:
            raise ValueError(
                "Lower limit should be smaller than upper limit")

        if not _CREATED[0]:
            from .. import _compat
            _compat._accum_fold[0] = _fold_coefficient
            _CREATED[0] = True
        return Expr.__new__(cls, min, max)

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

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        Add, Mul, Pow, S, Min, Max, And, FiniteSet, oo, zoo = _names()
        if isinstance(other, Expr):
            if isinstance(other, AccumBounds):
                return AccumBounds(
                    Add(self.min, other.min),
                    Add(self.max, other.max))
            if other is S.Infinity and self.min is S.NegativeInfinity or \
                    other is S.NegativeInfinity and self.max is S.Infinity:
                return AccumBounds(-oo, oo)
            elif other.is_extended_real:
                if self.min is S.NegativeInfinity and self.max is S.Infinity:
                    return AccumBounds(-oo, oo)
                elif self.min is S.NegativeInfinity:
                    return AccumBounds(-oo, self.max + other)
                elif self.max is S.Infinity:
                    return AccumBounds(self.min + other, oo)
                else:
                    return AccumBounds(Add(self.min, other), Add(self.max, other))
            return Add(self, other, evaluate=False)
        return NotImplemented

    __radd__ = __add__

    def __neg__(self):
        return AccumBounds(-self.max, -self.min)

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        Add, Mul, Pow, S, Min, Max, And, FiniteSet, oo, zoo = _names()
        if isinstance(other, Expr):
            if isinstance(other, AccumBounds):
                return AccumBounds(
                    Add(self.min, -other.max),
                    Add(self.max, -other.min))
            if other is S.NegativeInfinity and self.min is S.NegativeInfinity or \
                    other is S.Infinity and self.max is S.Infinity:
                return AccumBounds(-oo, oo)
            elif other.is_extended_real:
                if self.min is S.NegativeInfinity and self.max is S.Infinity:
                    return AccumBounds(-oo, oo)
                elif self.min is S.NegativeInfinity:
                    return AccumBounds(-oo, self.max - other)
                elif self.max is S.Infinity:
                    return AccumBounds(self.min - other, oo)
                else:
                    return AccumBounds(
                        Add(self.min, -other),
                        Add(self.max, -other))
            return Add(self, -other)
        return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __rsub__(self, other):
        return self.__neg__() + other

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        Add, Mul, Pow, S, Min, Max, And, FiniteSet, oo, zoo = _names()
        if self.args == (-oo, oo):
            return self
        if isinstance(other, Expr):
            if isinstance(other, AccumBounds):
                if other.args == (-oo, oo):
                    return other
                v = set()
                for a in self.args:
                    vi = other*a
                    v.update(vi.args or (vi,))
                return AccumBounds(Min(*v), Max(*v))
            if other is S.Infinity:
                if self.min.is_zero:
                    return AccumBounds(0, oo)
                if self.max.is_zero:
                    return AccumBounds(-oo, 0)
            if other is S.NegativeInfinity:
                if self.min.is_zero:
                    return AccumBounds(-oo, 0)
                if self.max.is_zero:
                    return AccumBounds(0, oo)
            if other.is_extended_real:
                if other.is_zero:
                    if self.max is S.Infinity:
                        return AccumBounds(0, oo)
                    if self.min is S.NegativeInfinity:
                        return AccumBounds(-oo, 0)
                    return S.Zero
                if other.is_extended_positive:
                    return AccumBounds(
                        Mul(self.min, other),
                        Mul(self.max, other))
                elif other.is_extended_negative:
                    return AccumBounds(
                        Mul(self.max, other),
                        Mul(self.min, other))
            if other.is_Order:
                return other
            return Mul(self, other, evaluate=False)
        return NotImplemented

    __rmul__ = __mul__

    @_sympifyit('other', NotImplemented)
    def __truediv__(self, other):
        Add, Mul, Pow, S, Min, Max, And, FiniteSet, oo, zoo = _names()
        if isinstance(other, Expr):
            if isinstance(other, AccumBounds):
                if other.min.is_positive or other.max.is_negative:
                    return self * AccumBounds(1/other.max, 1/other.min)

                if (self.min.is_extended_nonpositive and self.max.is_extended_nonnegative and
                    other.min.is_extended_nonpositive and other.max.is_extended_nonnegative):
                    if self.min.is_zero and other.min.is_zero:
                        return AccumBounds(0, oo)
                    if self.max.is_zero and other.min.is_zero:
                        return AccumBounds(-oo, 0)
                    return AccumBounds(-oo, oo)

                if self.max.is_extended_negative:
                    if other.min.is_extended_negative:
                        if other.max.is_zero:
                            return AccumBounds(self.max / other.min, oo)
                        if other.max.is_extended_positive:
                            # Union(Interval(-oo, self.max/other.max),
                            #       Interval(self.max/other.min, oo))
                            return AccumBounds(-oo, oo)

                    if other.min.is_zero and other.max.is_extended_positive:
                        return AccumBounds(-oo, self.max / other.max)

                if self.min.is_extended_positive:
                    if other.min.is_extended_negative:
                        if other.max.is_zero:
                            return AccumBounds(-oo, self.min / other.min)
                        if other.max.is_extended_positive:
                            # Union(Interval(-oo, self.min/other.min),
                            #       Interval(self.min/other.max, oo))
                            return AccumBounds(-oo, oo)

                    if other.min.is_zero and other.max.is_extended_positive:
                        return AccumBounds(self.min / other.max, oo)

            elif other.is_extended_real:
                if other in (S.Infinity, S.NegativeInfinity):
                    if self == AccumBounds(-oo, oo):
                        return AccumBounds(-oo, oo)
                    if self.max is S.Infinity:
                        return AccumBounds(Min(0, other), Max(0, other))
                    if self.min is S.NegativeInfinity:
                        return AccumBounds(Min(0, -other), Max(0, -other))
                if other.is_extended_positive:
                    return AccumBounds(self.min / other, self.max / other)
                elif other.is_extended_negative:
                    return AccumBounds(self.max / other, self.min / other)
            if (1 / other) is S.ComplexInfinity:
                return Mul(self, 1 / other, evaluate=False)
            else:
                return _flat_mul(self, 1 / other)

        return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __rtruediv__(self, other):
        Add, Mul, Pow, S, Min, Max, And, FiniteSet, oo, zoo = _names()
        if isinstance(other, Expr):
            if other.is_extended_real:
                if other.is_zero:
                    return S.Zero
                if (self.min.is_extended_nonpositive and self.max.is_extended_nonnegative):
                    if self.min.is_zero:
                        if other.is_extended_positive:
                            return AccumBounds(Mul(other, 1 / self.max), oo)
                        if other.is_extended_negative:
                            return AccumBounds(-oo, Mul(other, 1 / self.max))
                    if self.max.is_zero:
                        if other.is_extended_positive:
                            return AccumBounds(-oo, Mul(other, 1 / self.min))
                        if other.is_extended_negative:
                            return AccumBounds(Mul(other, 1 / self.min), oo)
                    return AccumBounds(-oo, oo)
                else:
                    return AccumBounds(Min(other / self.min, other / self.max),
                                       Max(other / self.min, other / self.max))
            return Mul(other, 1 / self, evaluate=False)
        else:
            return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __pow__(self, other):
        Add, Mul, Pow, S, Min, Max, And, FiniteSet, oo, zoo = _names()
        if isinstance(other, Expr):
            if other is S.Infinity:
                if self.min.is_extended_nonnegative:
                    if self.max < 1:
                        return S.Zero
                    if self.min > 1:
                        return S.Infinity
                    return AccumBounds(0, oo)
                elif self.max.is_extended_negative:
                    if self.min > -1:
                        return S.Zero
                    if self.max < -1:
                        return zoo
                    return S.NaN
                else:
                    if self.min > -1:
                        if self.max < 1:
                            return S.Zero
                        return AccumBounds(0, oo)
                    return AccumBounds(-oo, oo)

            if other is S.NegativeInfinity:
                return (1/self)**oo

            # generically true
            if (self.max - self.min).is_nonnegative:
                # well defined
                if self.min.is_nonnegative:
                    # no 0 to worry about
                    if other.is_nonnegative:
                        # no infinity to worry about
                        return self.func(self.min**other, self.max**other)

            if other.is_zero:
                return S.One  # x**0 = 1

            if other.is_Integer or other.is_integer:
                if self.min.is_extended_positive:
                    return AccumBounds(
                        Min(self.min**other, self.max**other),
                        Max(self.min**other, self.max**other))
                elif self.max.is_extended_negative:
                    return AccumBounds(
                        Min(self.max**other, self.min**other),
                        Max(self.max**other, self.min**other))

                if other % 2 == 0:
                    if other.is_extended_negative:
                        if self.min.is_zero:
                            return AccumBounds(self.max**other, oo)
                        if self.max.is_zero:
                            return AccumBounds(self.min**other, oo)
                        return (1/self)**(-other)
                    return AccumBounds(
                        S.Zero, Max(self.min**other, self.max**other))
                elif other % 2 == 1:
                    if other.is_extended_negative:
                        if self.min.is_zero:
                            return AccumBounds(self.max**other, oo)
                        if self.max.is_zero:
                            return AccumBounds(-oo, self.min**other)
                        return (1/self)**(-other)
                    return AccumBounds(self.min**other, self.max**other)

            # non-integer exponent
            # 0**neg or neg**frac yields complex
            if (other.is_number or other.is_rational) and (
                    self.min.is_extended_nonnegative or (
                    other.is_extended_nonnegative and
                    self.min.is_extended_nonnegative)):
                num, den = other.as_numer_denom()
                if num is S.One:
                    return AccumBounds(*[i**(1/den) for i in self.args])

                elif den is not S.One:  # e.g. if other is not Float
                    return (self**num)**(1/den)  # ok for non-negative base

            if isinstance(other, AccumBounds):
                if (self.min.is_extended_positive or
                        self.min.is_extended_nonnegative and
                        other.min.is_extended_nonnegative):
                    p = [self**i for i in other.args]
                    if not any(i.is_Pow for i in p):
                        a = [j for i in p for j in i.args or (i,)]
                        try:
                            return self.func(min(a), max(a))
                        except TypeError:  # can't sort
                            pass

            return Pow(self, other, evaluate=False)

        return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __rpow__(self, other):
        Add, Mul, Pow, S, Min, Max, And, FiniteSet, oo, zoo = _names()
        if other.is_real and other.is_extended_nonnegative and (
                self.max - self.min).is_extended_positive:
            if other is S.One:
                return S.One
            if other.is_extended_positive:
                a, b = [other**i for i in self.args]
                if min(a, b) != a:
                    a, b = b, a
                return self.func(a, b)
            if other.is_zero:
                if self.min.is_zero:
                    return self.func(0, 1)
                if self.min.is_extended_positive:
                    return S.Zero

        return Pow(other, self, evaluate=False)

    def __abs__(self):
        Add, Mul, Pow, S, Min, Max, And, FiniteSet, oo, zoo = _names()
        if self.max.is_extended_negative:
            return self.__neg__()
        elif self.min.is_extended_negative:
            return AccumBounds(S.Zero, Max(abs(self.min), self.max))
        else:
            return self

    def __contains__(self, other):
        """``True`` if ``other`` is contained in self, ``False`` if not,
        otherwise TypeError is raised."""
        Add, Mul, Pow, S, Min, Max, And, FiniteSet, oo, zoo = _names()
        other = _sympify(other)

        if other in (S.Infinity, S.NegativeInfinity):
            if self.min is S.NegativeInfinity or self.max is S.Infinity:
                return True
            return False

        rv = And(self.min <= other, self.max >= other)
        if rv not in (True, False):
            raise TypeError("input failed to evaluate")
        return rv

    def intersection(self, other):
        """The intersection of ``self`` and ``other`` (an AccumulationBounds
        or a FiniteSet)."""
        Add, Mul, Pow, S, Min, Max, And, FiniteSet, oo, zoo = _names()
        if not isinstance(other, (AccumBounds, FiniteSet)):
            raise TypeError(
                "Input must be AccumulationBounds or FiniteSet object")

        if isinstance(other, FiniteSet):
            fin_set = S.EmptySet
            for i in other:
                if i in self:
                    fin_set = fin_set + FiniteSet(i)
            return fin_set

        if self.max < other.min or self.min > other.max:
            return S.EmptySet

        if self.min <= other.min:
            if self.max <= other.max:
                return AccumBounds(other.min, self.max)
            if self.max > other.max:
                return other

        if other.min <= self.min:
            if other.max < self.max:
                return AccumBounds(self.min, other.max)
            if other.max > self.max:
                return self

    @property
    def contains(self):
        # The reference class has no ``contains`` (it is not a Set).
        raise AttributeError("'AccumulationBounds' object has no attribute 'contains'")

    def union(self, other):
        # The reference method (not a true union of disjoint bounds).
        Add, Mul, Pow, S, Min, Max, And, FiniteSet, oo, zoo = _names()
        if not isinstance(other, AccumBounds):
            raise TypeError(
                "Input must be AccumulationBounds or FiniteSet object")

        if self.min <= other.min and self.max >= other.min:
            return AccumBounds(self.min, Max(self.max, other.max))

        if other.min <= self.min and other.max >= self.min:
            return AccumBounds(other.min, Max(self.max, other.max))


def _is_ge(lhs, rhs):
    from ..core.relational import is_ge
    return is_ge(lhs, rhs)


def _is_lt(lhs, rhs):
    from ..core.relational import is_lt
    return is_lt(lhs, rhs)


def _is_le(lhs, rhs):
    from ..core.relational import is_le
    return is_le(lhs, rhs)


def _is_gt(lhs, rhs):
    from ..core.relational import is_gt
    return is_gt(lhs, rhs)


@dispatch(AccumulationBounds, AccumulationBounds)
def _eval_is_ge(lhs, rhs):  # noqa: F811
    if _is_ge(lhs.min, rhs.max):
        return True
    if _is_lt(lhs.max, rhs.min):
        return False


@dispatch(AccumulationBounds, Expr)  # noqa: F811
def _eval_is_ge(lhs, rhs):  # noqa: F811
    if not rhs.is_extended_real:
        raise TypeError(
            "Invalid comparison of %s %s" %
            (type(rhs), rhs))
    elif rhs.is_comparable:
        if _is_ge(lhs.min, rhs):
            return True
        if _is_lt(lhs.max, rhs):
            return False


@dispatch(Expr, AccumulationBounds)  # noqa: F811
def _eval_is_ge(lhs, rhs):  # noqa: F811
    if not lhs.is_extended_real:
        raise TypeError(
            "Invalid comparison of %s %s" %
            (type(lhs), lhs))
    elif lhs.is_comparable:
        if _is_le(rhs.max, lhs):
            return True
        if _is_gt(rhs.min, lhs):
            return False


# setting an alias for AccumulationBounds
AccumBounds = AccumulationBounds


def _fold_coefficient(op, args):
    """SymPy ``Add.flatten`` and ``Mul.flatten`` for terms that include
    bounds: the Numbers and the bounds fold into the coefficient through the
    bounds arithmetic, and the coefficient leads the arguments. ``None``
    leaves the construction to the native engine."""
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
    from .. import _core, sympify, S
    flat = []
    for arg in args:
        arg = sympify(arg)
        if getattr(arg, nested, False) and arg.is_commutative:
            flat.extend(arg.args)
        else:
            flat.append(arg)
    coeff = S.Zero if op == 'Add' else S.One
    rest = []
    for arg in flat:
        if op == 'Add' and arg.is_Mul and type(arg.args[0]) is AccumulationBounds:
            # A Number times bounds distributed over a sum.
            arg = _fold_coefficient('Mul', arg.args) or arg
        if type(arg) is AccumulationBounds:
            coeff = arg.__add__(coeff) if op == 'Add' else arg.__mul__(coeff)
        elif arg.is_Number:
            coeff = coeff + arg if op == 'Add' else coeff*arg
        else:
            rest.append(arg)
        if not (type(coeff) is AccumulationBounds or coeff.is_Number):
            return None
    if coeff is S.NaN:
        return coeff
    if not rest:
        return coeff
    return _core._node(op, (coeff, *rest), True)
__all__ = ['AccumBounds', 'AccumulationBounds']
