"""Admission for the reference ``_tsolve`` steps that the native solver
does not take (SymPy 1.14's ``sympy.solvers.solvers._tsolve`` for a single
equation in one symbol, after the native solver raised
NotImplementedError).

The steps (``_invert``, factoring, ``logcombine``, ``f(x)**e = rhs``, the
bivariate/Lambert block, the posify retry) and SymPy's ``checksol``
filtering of the candidates are native: ``symox::solve::tsolve``. Only the
admission, the recursion guard and the ``inverse()`` of host function
classes stay here. Returns None when no step applies, so that the caller
re-raises.
"""


_INVERSES = (
    ('acos', 'cos'), ('acosh', 'cosh'), ('acot', 'cot'), ('acoth', 'coth'),
    ('acsc', 'csc'), ('acsch', 'csch'), ('asec', 'sec'), ('asech', 'sech'),
    ('asin', 'sin'), ('asinh', 'sinh'), ('atan', 'tan'), ('atanh', 'tanh'),
    ('conjugate', 'conjugate'), ('cot', 'acot'), ('coth', 'acoth'),
    ('erf', 'erfinv'), ('erfc', 'erfcinv'), ('erfcinv', 'erfc'),
    ('erfinv', 'erf'), ('exp', 'log'), ('exp_polar', 'log'), ('log', 'exp'),
    ('sinh', 'asinh'), ('tan', 'atan'), ('tanh', 'atanh'))


def _inverse(f):
    """The reference ``f.inverse()`` of the one-argument functions that
    define it (``sin``, ``cos``, ``sec``, ``csc``, ``cosh`` and ``LambertW``
    define none)."""
    import symox
    inverse = getattr(f, 'inverse', None)
    if inverse is not None:
        try:
            return inverse()
        except (TypeError, NotImplementedError):
            return None
    for name, inv in _INVERSES:
        F = getattr(symox, name, None)
        if F is not None and isinstance(F, type) and isinstance(f, F):
            return getattr(symox, inv)
    return None


def _invert(eq, sym):
    """SymPy's ``_invert`` for one symbol (``integer_power=False``).
    Returns ``(rhs, lhs)``; native (``symox::solve::tsolve::invert``)."""
    from .. import _core, sympify
    return _core._tsolve_invert(sympify(eq), sympify(sym))


_SAW = []


def tsolve_step(solve, f, symbols, flags):
    """Like the reference ``tsolve_saw`` flag, an equation already being
    handled here is not handled again (that would recurse forever)."""
    key = (f, tuple(symbols))
    if key in _SAW:
        return None
    _SAW.append(key)
    try:
        return _tsolve_step(solve, f, symbols, flags)
    finally:
        _SAW.pop()


def _tsolve_step(solve, f, symbols, flags):
    from .. import Equality, Expr, _core, sympify
    if isinstance(f, (list, tuple, set)) or len(symbols) != 1:
        return None
    if flags.get('dict') or flags.get('set'):
        return None
    sym = sympify(symbols[0])
    if not getattr(sym, 'is_Symbol', False):
        return None
    f = sympify(f)
    if isinstance(f, Equality):
        f = f.lhs - f.rhs
    if not isinstance(f, Expr) or not f.has(sym):
        return None
    # the steps, the recursive solves they make and SymPy's check of the
    # candidates run natively (``symox::solve::tsolve``)
    return _core._tsolve_step(solve, f, sym, dict(flags))
