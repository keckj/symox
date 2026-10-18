"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox._public_tools', 'check_assumptions failing_assumptions')

# Expose native assumption properties through the standard introspection API.
_assume_defined = frozenset(('algebraic antihermitian commutative complex composite '
    'even extended_negative extended_nonnegative extended_nonpositive '
    'extended_nonzero extended_positive extended_real finite hermitian imaginary '
    'infinite integer irrational negative noninteger nonnegative nonpositive '
    'nonzero odd polar positive prime rational real transcendental zero').split())


def assumptions(expr, _check=None):
    from symox import sympify
    expr = sympify(expr)
    if expr.is_Symbol:
        facts = expr.assumptions0
        return facts if _check is None else {key:facts[key] for key in _check if key in facts}
    names = _assume_defined if _check is None else _check
    result = {}
    for name in names:
        value = getattr(expr, 'is_'+name, None)
        if value is not None:result[name] = value
    return result


def common_assumptions(exprs, check=None):
    facts = [assumptions(expr, _check=check) for expr in exprs]
    if not facts:return {}
    return {name:value for name,value in facts[0].items()
            if all(name in other and other[name] == value for other in facts[1:])}

__all__ += ['assumptions', 'common_assumptions']
