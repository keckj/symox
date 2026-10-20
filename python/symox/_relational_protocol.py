"""SymPy's comparison protocol for operands the native kernel cannot decide.

The native kernel decides orderings between native expressions. Instances of
user classes (Torch's ``int_oo``, ``SingletonInt``) carry their facts,
numeric values and comparison rules in Python: ``_eval_is_ge``
registrations, ``is_comparable`` and ``_as_mpf_val``. When a relation over
such operands stays unevaluated natively, it is decided here exactly as
SymPy's ``Relational.__new__`` and ``is_ge`` do (sympy/core/relational.py).
"""
from . import _core
from .multipledispatch import dispatch as _dispatch
from ._expr_facade import Expr as _Expr


@_dispatch(_Expr, _Expr)
def _eval_is_ge(lhs, rhs):
    return None


def is_comparable(self):
    """Return True if self can be computed to a real number
    (or already is a real number) with precision, else False."""
    if self.is_extended_real is False:
        return False
    if not self.is_number:
        return False
    # Numbers already evaluated are not re-evaluated: that would create
    # spurious precision.
    n, i = [p.evalf(2) if not p.is_Number else p for p in self.as_real_imag()]
    if not (i.is_Number and n.is_Number):
        return False
    if i:
        # Numbers with an imaginary part cannot be compared; with a
        # precision of 1 the question cannot be decided.
        return False
    return getattr(n, '_prec', -1) != 1


def _is_host(value):
    return isinstance(value, _core.Expr) and _core._host_contains(value)[1]


def _host_number_evalf(value, digits):
    # SymPy's evalf of a Number without an evalf table entry:
    # Float._new(self._as_mpf_val(prec), prec), at the working precision of
    # evalf (four guard bits) and rounded to the requested one.
    from mpmath import libmp, mpf
    from . import S, Float
    prec = libmp.dps_to_prec(digits)
    evaluate = getattr(type(value), '_eval_evalf', None)
    if evaluate is not None:
        return evaluate(value, prec + 4)
    as_mpf = getattr(value, '_as_mpf_val', None)
    if as_mpf is None:
        return None
    mpf_value = as_mpf(prec + 4)
    if mpf_value == libmp.finf:
        return S.Infinity
    if mpf_value == libmp.fninf:
        return S.NegativeInfinity
    if mpf_value == libmp.fnan:
        return S.NaN
    return Float(mpf(mpf_value), precision=prec)


def _evalf2(value):
    if _is_host(value) and value.is_Number:
        result = _host_number_evalf(value, 2)
        return value if result is None else result
    return value.evalf(2)


def _n2(a, b):
    """Return (a - b).evalf(2) if a and b are comparable, else None."""
    if a.is_comparable and b.is_comparable:
        dif = _evalf2(a - b)
        if dif.is_comparable:
            return dif


def is_ge(lhs, rhs):
    """SymPy's ``is_ge`` without assumptions."""
    from . import S
    if not (isinstance(lhs, _Expr) and isinstance(rhs, _Expr)):
        raise TypeError("Can only compare inequalities with Expr")
    retval = _eval_is_ge(lhs, rhs)
    if retval is not None:
        return retval
    n2 = _n2(lhs, rhs)
    if n2 is not None:
        # Float comparison for infinity avoids an infinite recursion.
        if n2 in (S.Infinity, S.NegativeInfinity):
            n2 = float(n2)
        return n2 >= 0
    if lhs.is_extended_real and rhs.is_extended_real:
        if (lhs.is_infinite and lhs.is_extended_positive) or (rhs.is_infinite and rhs.is_extended_negative):
            return True
        diff = lhs - rhs
        if diff is not S.NaN:
            rv = diff.is_extended_nonnegative
            if rv is not None:
                return rv
    return None


def _fuzzy_not(value):
    return None if value is None else not value


def host_inequality(op, lhs, rhs):
    """Truth value of the ordering ``lhs op rhs``, or None when undecided.

    Called by the native kernel for operands holding user class instances,
    after its own evaluation left the relation unevaluated.
    """
    from . import S
    for me in (lhs, rhs):
        if me.is_extended_real is False:
            raise TypeError("Invalid comparison of non-real %s" % me)
        if me is S.NaN:
            raise TypeError("Invalid NaN comparison")
    if op == '>=':
        return is_ge(lhs, rhs)
    if op == '<=':
        return is_ge(rhs, lhs)
    if op == '>':
        return _fuzzy_not(is_ge(rhs, lhs))
    return _fuzzy_not(is_ge(lhs, rhs))


def install():
    _core.Expr.is_comparable = property(is_comparable)
