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


def _host_free(*values):
    """No host (user) class anywhere in ``values``: the reference
    constructors could defer to a user ``_eval_*`` hook."""
    from ... import _host
    if not _host.host_classes_defined():
        return True
    try:
        return not any(_core._host_contains(v)[0] for v in values)
    except Exception:
        return False


def _applied_function(func, x):
    return getattr(func, 'args', None) == (x,) and getattr(func, 'is_Function', False)


def undetermined_trialset(expr, x, func, eq_homogeneous, zero):
    """``_undetermined_coefficients_match`` when its test holds and the
    trial set is native; ``eq_homogeneous is zero`` selects the reference's
    branch without multiplication by powers of ``x``."""
    kernel = getattr(_core, 'ode_uc_trialset', None)
    if kernel is None or not getattr(x, 'is_Symbol', False):
        return None
    shift = eq_homogeneous is not zero
    if shift and not _applied_function(func, x):
        return None
    if not _host_free(expr, eq_homogeneous):
        return None
    try:
        fx = func if shift else x
        terms = kernel(expr, x, fx, eq_homogeneous if shift else None)
    except Exception:
        return None
    if terms is None:
        return None
    return {'test': True, 'trialset': set(terms)}


def undetermined_solution(eq, func, order, match, trialset):
    """The right-hand side of ``_solve_undetermined_coefficients``."""
    kernel = getattr(_core, 'ode_uc_solve', None)
    x = _variable(func)
    if kernel is None or x is None or not isinstance(match, dict):
        return None
    try:
        gensols = match['list']
        gsol_rhs = match['sol'].rhs
        n = len(gensols)
        terms = list(trialset)
    except Exception:
        return None
    if not _host_free(eq, gsol_rhs, *terms):
        return None
    try:
        return kernel(eq, func, x, int(order), n, gsol_rhs, terms)
    except Exception:
        return None


def variation_of_parameters(func, roots, homogen_sol, order, match_obj, simplify_flag,
                            simplify_wronskian, simplify_psol):
    """The right-hand side of ``_solve_variation_of_parameters``; the two
    simplifications are the reference's, run as callbacks."""
    kernel = getattr(_core, 'ode_vop_solve', None)
    x = _variable(func)
    if kernel is None or x is None or not isinstance(match_obj, dict):
        return None
    try:
        forcing = match_obj[-1]
        lead = match_obj[order]
        rhs = homogen_sol.rhs
        roots = list(roots)
        order = int(order)
    except Exception:
        return None
    if not _host_free(forcing, lead, rhs, *roots):
        return None
    called = []

    def wrap(f):
        def run(e):
            called.append(True)
            return f(e)
        return run
    try:
        return kernel(roots, x, forcing, lead, order, rhs, bool(simplify_flag),
                      wrap(simplify_wronskian), wrap(simplify_psol))
    except Exception:
        # an exception of the reference's simplifications is the reference's
        if called:
            raise
        return None
