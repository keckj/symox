"""Admission of the native characteristic-equation solutions of linear ODEs.

The characteristic polynomials, their roots, the fundamental solutions and
the collected general solutions are computed in Rust
(``symox::ode::constcoeff``). This module only checks the shape of the
inputs and converts the results; ``None`` means the kernel declined and the
reference body in ``nonhomogeneous`` must run unchanged.
"""

from ... import _core


def _integer_items(r):
    items = []
    for i in r.keys():
        if isinstance(i, str):
            continue
        if not isinstance(i, int) or isinstance(i, bool):
            return None
        if i >= 0:
            items.append((i, r[i]))
    return items


def _variable(func):
    args = getattr(func, 'args', ())
    if len(args) != 1 or not getattr(args[0], 'is_Symbol', False):
        return None
    return args[0]


def const_characteristic_sols(r, func, order):
    """``(gensols, collectterms)`` of ``_get_const_characteristic_eq_sols``."""
    kernel = getattr(_core, 'ode_const_char_sols', None)
    x = _variable(func)
    if kernel is None or x is None or not isinstance(r, dict):
        return None
    items = _integer_items(r)
    if items is None:
        return None
    try:
        native = kernel(items, x, int(order))
    except Exception:
        return None
    if native is None:
        return None
    gensols, collectterms = native
    return list(gensols), list(collectterms)


def simplified_sol(sol, func, collectterms):
    """The simplified right-hand side of ``_get_simplified_sol``."""
    kernel = getattr(_core, 'ode_const_simplified_sol', None)
    x = _variable(func)
    if kernel is None or x is None or len(sol) != 1:
        return None
    if getattr(sol[0], 'lhs', None) != func:
        return None
    try:
        return kernel(sol[0].rhs, x, list(collectterms))
    except Exception:
        return None


def euler_characteristic_sols(eq, func, match_obj, numbered_constants):
    """``(gsol_rhs, gensols)`` of ``_get_euler_characteristic_eq_sols``."""
    kernel = getattr(_core, 'ode_euler_char_sols', None)
    x = _variable(func)
    if kernel is None or x is None or not isinstance(match_obj, dict):
        return None
    items = _integer_items(match_obj)
    if not items:
        return None
    try:
        order = max(i for i, _ in items)
        constants = list(numbered_constants(eq, num=order*2))
        native = kernel(items, x, constants)
    except Exception:
        return None
    if native is None:
        return None
    rhs, gensols = native
    return rhs, list(gensols)
