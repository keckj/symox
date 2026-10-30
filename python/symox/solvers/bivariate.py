"""SymPy 1.14's ``sympy.solvers.bivariate``.

``_solve_lambert`` solves equations that become ``x = W(...)`` forms after
taking logarithms (``x*exp(x**2) = 1``, ``log(x) + x**2 = 0``,
``x**a*exp(x) = 1``) and ``bivariate_type`` recognizes ``f(x*y)``,
``f(a*x + b*y)`` and ``f(a*x*y + b*y)``. They are used by the ``_tsolve``
fallback of ``solve`` after the native solver raised NotImplementedError.

All of them are native (``symox::solve::tsolve``); these are wrappers.
"""


def _core():
    from .. import _core
    return _core


def _s(e):
    from .. import sympify
    return sympify(e)


def _mexpand(expr, recursive=False):
    if expr is None:
        return None
    return _core()._bivariate_mexpand(_s(expr), recursive)


def _filtered_gens(poly, symbol):
    """Generators of ``poly`` that contain ``symbol``; of ``g`` and ``1/g``
    only one is kept."""
    return set(_core()._bivariate_filtered_gens(list(poly.gens), _s(symbol)))


def _mostfunc(lhs, func, X=None):
    """The term of ``lhs`` with ``func`` that holds the most ``func``."""
    from .. import Pow, exp, log
    name = {log: 'log', exp: 'exp', Pow: 'Pow'}.get(func)
    if name is None:
        raise ValueError('unsupported class %s' % (func,))
    return _core()._bivariate_mostfunc(
        _s(lhs), name, None if not X else _s(X))


def _linab(arg, symbol):
    """``(a, b, x)`` with ``arg = a*x + b``, ``x`` the symbol-dependent part."""
    return _core()._bivariate_linab(_s(arg), _s(symbol))


def _lambert(eq, x):
    """Solutions of ``a*log(b*X + c) + d*X + f = 0`` with ``X = g(x)``."""
    return _core()._bivariate_lambert(_s(eq), _s(x))


def _solve_lambert(f, symbol, gens):
    """SymPy's ``_solve_lambert``: raises NotImplementedError when ``f``
    has no Lambert form."""
    return _core()._bivariate_solve_lambert(
        _s(f), _s(symbol), [_s(g) for g in gens])


def bivariate_type(f, x, y, *, first=True):
    """``(g, p, u)`` such that ``f(x, y) = p(u)`` with ``u = g(x, y)`` one
    of ``x*y``, ``a*x + b*y`` or ``a*x*y + b*y``, else None."""
    if not first:
        f = f.as_expr()
    return _core()._bivariate_type(_s(f), _s(x), _s(y))
