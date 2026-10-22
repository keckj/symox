"""The ``frac``, ``numer``, ``denom`` and ``modulus`` hints of ``expand``.

Follows SymPy 1.14's ``Expr.expand`` (``sympy/core/expr.py``, BSD-3-Clause,
see LICENSE-SYMPY): the fraction hints expand the parts returned by
``fraction`` and ``modulus`` reduces the rational coefficients of the
expanded terms.
"""
import sys
from functools import wraps

from . import _core

_FRACTION_HINTS = ('frac', 'numer', 'denom')


def _wanted(hints):
    return any(name in hints for name in _FRACTION_HINTS + ('modulus', 'exact'))


def _expand(expand, expr, deep, hints):
    from . import Add, sympify
    from .simplify.radsimp import fraction
    hints = dict(hints)
    modulus = hints.pop('modulus', None)
    exact = hints.pop('exact', False)
    frac = hints.pop('frac', False)
    denom = hints.pop('denom', False)
    numer = hints.pop('numer', False)

    def again(value):
        return _expand_with(expand, value, deep, modulus, exact, hints)

    if frac:
        n, d = [again(a) for a in fraction(expr, exact)]
        return n/d
    if denom:
        n, d = fraction(expr, exact)
        return n/again(d)
    if numer:
        n, d = fraction(expr, exact)
        return again(n)/d
    return _expand_with(expand, expr, deep, modulus, exact, hints)


def _expand_with(expand, expr, deep, modulus, exact, hints):
    from . import Add, sympify
    expr = expand(expr, deep, hints)
    if modulus is not None:
        modulus = sympify(modulus)
        if not modulus.is_Integer or modulus <= 0:
            raise ValueError("modulus must be a positive integer, got %s" % modulus)
        terms = []
        for term in Add.make_args(expr):
            coeff, tail = term.as_coeff_Mul(rational=True)
            coeff %= modulus
            if coeff:
                terms.append(coeff*tail)
        expr = Add(*terms)
    return expr


def install(namespace):
    original_root = namespace['expand']
    original_method = _core.Expr.expand

    def method_expand(value, deep, hints):
        return original_method(value, deep=deep, **hints)

    @wraps(original_method)
    def method(self, deep=True, **hints):
        if hints and _wanted(hints):
            return _expand(method_expand, self, deep, hints)
        return original_method(self, deep=deep, **hints)

    @wraps(original_root)
    def expand(e, deep=True, **hints):
        if hints and _wanted(hints):
            from . import sympify
            value = sympify(e)
            if isinstance(value, _core.Expr) and not getattr(value, 'is_Matrix', False):
                return _expand(method_expand, value, deep, hints)
        return original_root(e, deep=deep, **hints)

    namespace['expand'] = expand
    _core.Expr.expand = method
    for name, module in tuple(sys.modules.items()):
        if (name.startswith('symox.') and name != 'symox._core' and module is not None
                and module.__dict__.get('expand') is original_root):
            module.expand = expand
