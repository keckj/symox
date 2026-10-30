"""Reference ``solveset`` steps that the native solver does not take.

``modular_step`` admits ``A - Mod(B, C) = 0`` to the native port of SymPy
1.14's ``_solve_modular``/``_invert_modular``; ``radical_step`` (SymPy's
``_solve_radical`` after its ``unrad`` test, every candidate checked with
the native ``checksol``, an undecided candidate kept in a ConditionSet),
``negative_base_step`` and ``negative_base_guard`` run natively
(``crates/symox-core/src/sets/solveset_radical.rs``). The recursive
``solveset``, ``imageset`` and ``intersect`` calls are Python callbacks.
"""


def modular_step(f, symbol, domain, solver=None):
    """The reference result for ``A - Mod(B, C) = 0``, or None when ``f`` is
    not modular in ``symbol``. ``_solve_modular`` and ``_invert_modular``
    run natively (``crates/symox-core/src/sets/solveset_modular.rs``) and
    raise the reference TypeError, NotInvertible or AttributeError."""
    from .. import Expr, _core, _unwrap, _wrap_set, sympify
    f = sympify(f)
    if not isinstance(f, Expr) and not getattr(f, 'is_Relational', False):
        return None
    result = _core._solveset_modular(f, symbol, _unwrap(domain))
    return None if result is None else _wrap_set(result)


def _admit(f, symbol, domain):
    """The equation as an Expr ``f`` with a Symbol and a Set domain, or None."""
    from .. import Equality, Expr, Set, sympify
    if not isinstance(domain, Set) or not getattr(symbol, 'is_Symbol', False):
        return None
    f = sympify(f)
    if isinstance(f, Equality):
        f = f.lhs - f.rhs
    return f if isinstance(f, Expr) else None


def _callbacks(solver):
    from .. import _unwrap, _wrap_set, imageset, solve
    return (lambda f, x, d: _unwrap(solver(f, x, _wrap_set(d))),
            lambda lam, base: _unwrap(imageset(lam, _wrap_set(base))),
            lambda a, b: _unwrap(_wrap_set(a).intersect(_wrap_set(b))),
            solve)


def _unrad_fallback(f, x):
    """The established resultant elimination where SymPy's ``unrad`` raises
    NotImplementedError (``_unrad_resultant``), with core expressions."""
    from .. import _unwrap
    from ._algebraic import _unrad_resultant
    r = _unrad_resultant(f, x)
    return None if r is None else (_unwrap(r[0]), [_unwrap(c) for c in r[1]])


def radical_step(f, symbol, domain, solver):
    """SymPy's ``_solve_radical`` for an equation that the native solver
    leaves as a ConditionSet (``radical_step`` in
    ``crates/symox-core/src/sets/solveset_radical.rs``). Returns None to
    keep the native result."""
    from .. import _core, _unwrap, _wrap_set
    f = _admit(f, symbol, domain)
    if f is None:
        return None
    r = _core._solveset_radical_step(f, symbol, _unwrap(domain), *_callbacks(solver),
                                     _unrad_fallback)
    return None if r is None else _wrap_set(r)


def negative_base_step(f, symbol, domain, solver):
    """``_invert_real`` of an equation whose symbol sits in the exponent of a
    power with a negative constant base (``(-2)**x = 4``) on a real domain
    (native ``negative_base_step``). Returns None to keep the native
    result."""
    from .. import _core, _unwrap, _wrap_set
    f = _admit(f, symbol, domain)
    if f is None:
        return None
    r = _core._solveset_negative_base_step(f, symbol, _unwrap(domain), *_callbacks(solver))
    return None if r is None else _wrap_set(r)


def negative_base_guard(f, symbol, domain, result):
    """The native result becomes a ConditionSet when the principal logarithm
    may have missed the integer exponents of a negative base on a real
    domain (native ``negative_base_guard``)."""
    from .. import ConditionSet, _core, _unwrap, _wrap_set
    if isinstance(result, ConditionSet):
        return result
    f0 = _admit(f, symbol, domain)
    if f0 is None:
        return result
    r = _wrap_set(_core._solveset_negative_base_guard(f0, symbol, _unwrap(domain), _unwrap(result)))
    return result if r == result else r
