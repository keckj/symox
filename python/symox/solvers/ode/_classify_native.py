"""Admission of the native first-order ``classify_ode`` matches.

The decompositions are computed in Rust (``symox::ode::classify``); this
module only decides which solvers may use them and maps the native values
onto each solver's own wilds. A ``None`` entry from the kernel means the
reference matcher must run.
"""

from ... import _core

_MISSING = object()
_ROLES = None


def _roles():
    global _ROLES
    if _ROLES is None:
        from . import single as s
        _ROLES = {
            s.Separable._equation: ('de', (0, 1)),
            s.FirstExact._equation: ('de', (0, 1)),
            s.AlmostLinear._equation: ('de', (1, 0)),
            s.HomogeneousCoeffSubsDepDivIndep._equation: ('de', (0, 1)),
            s.HomogeneousCoeffSubsIndepDivDep._equation: ('de', (0, 1)),
            s.LinearCoefficients._equation: ('de', (0, 1)),
            s.FirstLinear._equation: ('linear', (0, 1)),
            s.Bernoulli._equation: ('bernoulli', (0, 1, 2)),
            s.RiccatiSpecial._equation: ('riccati_special', ()),
        }
        _ROLES = {k: v + (_original_wilds(s, k),) for k, v in _ROLES.items()}
    return _ROLES


def _original_wilds(s, equation):
    for cls in (s.Separable, s.FirstExact, s.AlmostLinear, s.HomogeneousCoeffSubsDepDivIndep,
                s.HomogeneousCoeffSubsIndepDivDep, s.LinearCoefficients, s.FirstLinear,
                s.Bernoulli, s.RiccatiSpecial):
        if cls._equation is equation:
            return cls._wilds
    return None


def _first_order(problem):
    cached = problem.__dict__.get('_native_first_order', _MISSING)
    if cached is not _MISSING:
        return cached
    result = None
    native = getattr(_core, 'ode_first_order_matches', None)
    if native is not None and problem.order == 1:
        eq = problem.eq_expanded
        if getattr(eq, 'is_commutative', False) and isinstance(eq, _core.Expr):
            from symox.polys.polytools import cancel
            from symox.core.function import expand
            from symox._operation_counts import count_ops
            fx = problem.func
            x = problem.sym
            try:
                result = native(eq, fx, x, fx.diff(x), cancel, expand, count_ops)
            except Exception:
                result = None
    problem.__dict__['_native_first_order'] = result
    return result


def solver_match(solver):
    """The native match of ``solver``: a dict over its wilds, ``False`` for no
    match, ``None`` when the reference matcher must run."""
    cls = type(solver)
    role = _roles().get(cls._equation)
    if role is None or cls._wilds is not role[2]:
        return None
    problem = solver.ode_problem
    if problem.order != 1:
        return None
    native = _first_order(problem)
    if native is None:
        return None
    value = native.get(role[0])
    if value is None or value is False:
        return value
    if role[0] == 'riccati_special':
        return None
    wilds = solver.wilds()
    return {wilds[i]: value[j] for i, j in enumerate(role[1]) if j < len(value)}


def split_de(collected, df, d, e, df2):
    """``collected.match(d + e*df)`` for the order-one power series match:
    a dict, ``None`` for no match, ``_MISSING`` when not decided."""
    native = getattr(_core, 'ode_split_de', None)
    if native is None or not isinstance(collected, _core.Expr) or not getattr(collected, 'is_commutative', False):
        return _MISSING
    from symox._operation_counts import count_ops
    try:
        value = native(collected, df, [df2], count_ops)
    except Exception:
        return _MISSING
    if value is None:
        return _MISSING
    if value is False:
        return None
    return {d: value[0], e: value[1]}


def ordered_hints(matching_hints):
    """``[i for i in allhints if i in matching_hints]``."""
    native = getattr(_core, 'ode_ordered_hints', None)
    keys = [k for k in matching_hints if isinstance(k, str)]
    if native is None:
        return None
    return list(native(keys))


def riccati_match(eq, f, x):
    """``match_riccati(eq, f, x)`` decided natively: ``(True, [b0, b1, b2])``,
    ``(False, [])``, or ``None`` when the reference must run."""
    native = getattr(_core, 'ode_match_riccati', None)
    if native is None or not isinstance(eq, _core.Expr) or not getattr(eq, 'is_commutative', False):
        return None
    from symox.polys.polytools import cancel
    from symox.core.function import expand
    fx = f(x)
    try:
        value = native(eq, fx, x, fx.diff(x), cancel, expand)
    except Exception:
        return None
    if value is None:
        return None
    if value is False:
        return False, []
    return True, list(value)


def _linear_args(eq, func):
    if not isinstance(eq, _core.Expr) or not getattr(eq, 'is_commutative', False):
        return None
    if not isinstance(func, _core.Expr) or len(func.args) != 1:
        return None
    return func.args[0]


def linear_coefficients(eq, func, order):
    """``get_linear_coefficients(eq, func, order)`` decided natively: the
    reference dict, ``None`` (not linear), or ``_MISSING`` (not decided)."""
    native = getattr(_core, 'ode_linear_coefficients', None)
    x = _linear_args(eq, func)
    if native is None or x is None or not isinstance(order, int) or order < 0:
        return _MISSING
    try:
        value = native(eq, func, x, order)
    except Exception:
        return _MISSING
    if value is None:
        return _MISSING
    if value is False:
        return None
    terms = {i: value[i + 1] for i in range(order + 1)}
    terms[-1] = value[0]
    return terms


def nth_linear_match(solver, eq):
    """The ``nth_linear_*`` matcher of ``solver`` on ``eq`` decided natively:
    ``(True, r)`` with the solver's ``self.r``, ``(False, None)``, or ``None``
    when the reference matcher must run."""
    native = getattr(_core, 'ode_nth_linear_match', None)
    problem = solver.ode_problem
    func = problem.func
    order = problem.order
    x = _linear_args(eq, func)
    if native is None or x is None or x != problem.sym or not isinstance(order, int) or order < 0:
        return None
    try:
        value = native(eq, func, x, order, solver.hint)
    except Exception:
        return None
    if value is None:
        return None
    if value is False:
        return False, None
    r = {i: value[i + 1] for i in range(order + 1)}
    r[-1] = value[0]
    return True, r


class LazyTrialset:
    """``self.trialset`` of the constant-coefficient undetermined-coefficients
    solver after a native match: the reference
    ``_undetermined_coefficients_match(...)['trialset']``, computed on first
    use by the solution step. An assignment (the reference matcher) stores
    the value as before."""

    def __set_name__(self, owner, name):
        self.name = '_' + name

    def __get__(self, obj, owner=None):
        if obj is None:
            return self
        d = obj.__dict__
        if self.name not in d:
            args = d.get('_native_trial_args')
            if args is None:
                raise AttributeError('trialset')
            from .nonhomogeneous import _undetermined_coefficients_match
            d[self.name] = _undetermined_coefficients_match(*args)['trialset']
        return d[self.name]

    def __set__(self, obj, value):
        obj.__dict__[self.name] = value

