"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'Eq Equality Ge GreaterThan Gt Le LessThan Lt Ne Rel Relational StrictGreaterThan StrictLessThan Unequality')


def is_ge(lhs, rhs, assumptions=None):
    """Return a three-valued extended-real comparison without coercing inputs."""
    from symox import Expr, Ge, Q, ask, true, false, nan
    if not isinstance(lhs, Expr) or not isinstance(rhs, Expr):
        raise TypeError('Can only compare inequalities with Expr')
    if lhs == nan or rhs == nan:
        return None
    def fact(value, name):
        if assumptions is None:
            return getattr(value, 'is_'+name, None)
        return ask(getattr(Q, name)(value), assumptions)
    if fact(lhs, 'extended_real') is not True or fact(rhs, 'extended_real') is not True:
        return None
    result = Ge(lhs, rhs)
    # The numeric branch returns a symbolic Boolean in SymPy; assumption
    # deductions return Python bools. Preserve that public distinction.
    def comparable(value):
        return bool(value.is_number and value.is_extended_real and value.is_finite)
    if (comparable(lhs) and comparable(rhs) and comparable(lhs-rhs)
            and (result == true or result == false)):
        return result
    if result == true:
        return True
    if result == false:
        return False
    return fact(lhs-rhs, 'extended_nonnegative')


def is_le(lhs, rhs, assumptions=None):
    return is_ge(rhs, lhs, assumptions)


def is_lt(lhs, rhs, assumptions=None):
    from .logic import fuzzy_not
    return fuzzy_not(is_ge(lhs, rhs, assumptions))


def is_gt(lhs, rhs, assumptions=None):
    from .logic import fuzzy_not
    return fuzzy_not(is_le(lhs, rhs, assumptions))


__all__ += ['is_ge', 'is_le', 'is_lt', 'is_gt']
