"""Exact Dirac expansion through the native polynomial and root engines."""
from functools import wraps
import sys

from . import _core


def _delta_node(node, wrt, max_terms):
    from . import (Dummy, Poly, DiracDelta, solve, sympify, S,
                   Add)
    argument = node.args[0]
    if wrt is None:
        variables = node.free_symbols
        if len(variables) != 1:
            raise TypeError("The 'wrt' keyword is required when expanding a DiracDelta with zero or several free symbols")
        variable = next(iter(variables))
    else:
        variable = sympify(wrt)
        if not getattr(variable, 'is_Symbol', False):
            raise ValueError('DiracDelta expansion requires a symbol as wrt')
    if not argument.has(variable) or len(node.args) > 1 and node.args[1] != 0:
        return node
    if _core._exact_assumption(variable, 'real') is False:
        return node
    # Floating roots require a separate precision contract. Never pass them
    # through the polynomial solver's approximate root path here.
    if any(value.is_Float for value in argument.atoms() if isinstance(value, _core.Expr)):
        return node
    if argument == variable:
        return node
    slope = argument.diff(variable)
    if not slope.has(variable) and _core._exact_assumption(slope, 'real') is True and _core._exact_assumption(slope, 'zero') is False:
        offset = argument.subs(variable, 0)
        root = -offset/slope
        if _core._exact_assumption(root, 'real') is True:
            denominator = _core._distribution_abs_at(slope, variable, 0)
            return node if denominator is None else DiracDelta(variable-root)/denominator
        return node
    witness = Dummy('delta_expansion')
    if argument.is_polynomial(variable) is not True:
        return node
    polynomial = Poly(argument.xreplace({variable: witness}), witness)
    if any(coefficient.is_Rational is not True for coefficient in polynomial.all_coeffs()):
        return node
    degree = polynomial.degree()
    if degree <= 0:
        return node
    derivative = polynomial.diff()
    # A nonconstant gcd detects multiple roots before any division by g'(r).
    if polynomial.gcd(derivative).degree() > 0:
        return node
    roots = solve(polynomial.as_expr(), witness, cubics=False, quartics=False)
    if not isinstance(roots, list) or len(roots) != degree:
        return node
    real = []
    for root in roots:
        status = _core._exact_assumption(root, 'real')
        if status is None:
            return node
        if status:
            real.append(root)
    if len(real) > max_terms:
        raise ValueError('DiracDelta expansion exceeds max_terms')
    terms = []
    for root in real:
        denominator = _core._distribution_abs_at(derivative.as_expr(), witness, root)
        if denominator is None:
            return node
        terms.append(DiracDelta(variable-root)/denominator)
    return Add(*terms) if terms else S.Zero


def _expand_deltas(expression, wrt, deep, max_terms):
    from . import MatrixBase
    if isinstance(expression, MatrixBase):
        return expression.applyfunc(lambda value: _expand_deltas(value, wrt, deep, max_terms))
    def is_delta(value):
        return (getattr(value, 'class_name', None) == 'DiracDelta'
                and not getattr(value, 'is_UndefinedFunction', False))
    def expand_node(value):
        from .polys.polyerrors import PolynomialError
        try:
            return _delta_node(value, wrt, max_terms)
        except (PolynomialError, NotImplementedError):
            return value
    if not deep:
        return expand_node(expression) if is_delta(expression) else expression
    return expression.replace(is_delta, expand_node)


def install(namespace):
    original_root = namespace['expand']
    original_method = _core.Expr.expand

    def options(hints):
        hints = dict(hints)
        enabled = hints.pop('diracdelta', False)
        wrt = hints.pop('wrt', None)
        return hints, enabled, wrt

    @wraps(original_root)
    def expand(e, deep=True, **hints):
        rest, enabled, wrt = options(hints)
        result = original_root(e, deep=deep, **rest)
        return _expand_deltas(result, wrt, deep, rest.get('max_terms', 100_000)) if enabled else result

    @wraps(original_method)
    def method(self, deep=True, **hints):
        rest, enabled, wrt = options(hints)
        result = original_method(self, deep=deep, **rest)
        return _expand_deltas(result, wrt, deep, rest.get('max_terms', 100_000)) if enabled else result

    namespace['expand'] = expand
    _core.Expr.expand = method
    for name, module in tuple(sys.modules.items()):
        if (name.startswith('symox.') and name != 'symox._core' and module is not None
                and module.__dict__.get('expand') is original_root):
            module.expand = expand
