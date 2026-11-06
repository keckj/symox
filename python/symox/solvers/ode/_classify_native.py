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
            s.RiccatiSpecial._equation: ('riccati_special', (0, 1, 2, 3)),
            s.Liouville._equation: ('liouville', (0, 1, 2)),
            s.SecondNonlinearAutonomousConserved._equation: ('conserved', (0,)),
        }
        _ROLES = {k: v + (_original_wilds(s, k),) for k, v in _ROLES.items()}
    return _ROLES


def _original_wilds(s, equation):
    for cls in (s.Separable, s.FirstExact, s.AlmostLinear, s.HomogeneousCoeffSubsDepDivIndep,
                s.HomogeneousCoeffSubsIndepDivDep, s.LinearCoefficients, s.FirstLinear,
                s.Bernoulli, s.RiccatiSpecial, s.Liouville, s.SecondNonlinearAutonomousConserved):
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


def _second_order(problem):
    cached = problem.__dict__.get('_native_second_order', _MISSING)
    if cached is not _MISSING:
        return cached
    result = None
    native = getattr(_core, 'ode_second_order_patterns', None)
    if native is not None and problem.order == 2:
        eq = problem.eq_expanded
        if getattr(eq, 'is_commutative', False) and isinstance(eq, _core.Expr):
            from symox.polys.polytools import cancel
            from symox.core.function import expand
            try:
                result = native(eq, problem.func, problem.sym, cancel, expand)
            except Exception:
                result = None
    problem.__dict__['_native_second_order'] = result
    return result


def solver_match(solver):
    """The native match of ``solver``: a dict over its wilds, ``False`` for no
    match, ``None`` when the reference matcher must run."""
    cls = type(solver)
    role = _roles().get(cls._equation)
    if role is None or cls._wilds is not role[2]:
        return None
    problem = solver.ode_problem
    if problem.order == 1:
        native = _first_order(problem)
    elif problem.order == 2 and role[0] in ('liouville', 'conserved'):
        native = _second_order(problem)
    else:
        return None
    if native is None:
        return None
    value = native.get(role[0])
    if value is None or value is False:
        return value
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


def _second_args(eq, func):
    x = _linear_args(eq, func)
    if x is None or not isinstance(x, _core.Expr):
        return None
    return x


def match_2nd_hypergeometric(eq, func):
    """``match_2nd_hypergeometric(eq, func)`` decided natively: ``[A, B]``,
    ``[]``, or ``_MISSING`` (not decided)."""
    from symox.polys.polytools import cancel
    native = getattr(_core, 'ode_match_2nd_hypergeometric', None)
    x = _second_args(eq, func)
    if native is None or x is None:
        return _MISSING
    try:
        value = native(eq, func, x, cancel)
    except Exception:
        return _MISSING
    if value is None:
        return _MISSING
    if value is False:
        return []
    return list(value)


def _hyper_hooks():
    from symox.polys.polytools import cancel, factor
    from symox.polys import roots
    from symox.core.function import expand
    from symox.simplify import simplify
    return cancel, factor, simplify, expand, roots


def hyper_2f1_match(A, B, func):
    """``equivalence_hypergeometric`` and ``match_2nd_2F1_hypergeometric``
    decided natively: the reference match object, ``None``, or ``_MISSING``
    (not decided)."""
    native = getattr(_core, 'ode_hyper_2f1_match', None)
    if native is None or not isinstance(func, _core.Expr) or len(func.args) != 1:
        return _MISSING
    if not (isinstance(A, _core.Expr) and isinstance(B, _core.Expr)):
        return _MISSING
    try:
        value = native(A, B, func.args[0], *_hyper_hooks())
    except Exception:
        return _MISSING
    if value is None:
        return _MISSING
    if value is False:
        return None
    a, b, c, k, mobius = value
    return {'a': a, 'b': b, 'c': c, 'k': k, 'mobius': mobius, 'type': "2F1"}


def hyper_2f1_solution(eq, func, match_object):
    """``get_sol_2F1_hypergeometric`` decided natively: ``Eq(f(x), rhs)``,
    ``None``, or ``_MISSING`` (not decided)."""
    native = getattr(_core, 'ode_hyper_2f1_solution', None)
    if native is None or not isinstance(func, _core.Expr) or len(func.args) != 1:
        return _MISSING
    params = tuple(match_object.get(key) for key in ('a', 'b', 'c', 'k', 'mobius'))
    A = match_object.get('A')
    if not all(isinstance(v, _core.Expr) for v in params + (A,)):
        return _MISSING
    from .ode import get_numbered_constants
    from symox.core.relational import Eq
    C0, C1 = get_numbered_constants(eq, num=2)
    try:
        value = native(params, A, func.args[0], C0, C1, *_hyper_hooks())
    except Exception:
        return _MISSING
    if value is None:
        return _MISSING
    if value is False:
        return None
    return Eq(func, value[0])


def airy_match(solver, eq):
    """The ``2nd_linear_airy`` matcher decided natively: ``(True, rn)``,
    ``(False, None)``, or ``None`` when the reference matcher must run."""
    from symox.polys.polytools import cancel
    native = getattr(_core, 'ode_airy_match', None)
    problem = solver.ode_problem
    order = problem.order
    x = _second_args(eq, problem.func)
    if native is None or x is None or x != problem.sym or not isinstance(order, int) or order < 0:
        return None
    try:
        value = native(eq, problem.func, x, order, cancel)
    except Exception:
        return None
    if value is None:
        return None
    if value is False:
        return False, None
    return True, {'b': value[0], 'm': value[1]}


def bessel_rejected(solver, eq):
    """``True`` when the ``2nd_linear_bessel`` matcher provably returns
    ``False``; otherwise the reference matcher must run."""
    native = getattr(_core, 'ode_bessel_rejected', None)
    problem = solver.ode_problem
    x = _second_args(eq, problem.func)
    if native is None or x is None or x != problem.sym:
        return False
    try:
        return native(eq, problem.func, x) is False
    except Exception:
        return False


def power_series2(reduced_eq, f, x, point, wilds, terms):
    """The order-two power series hints of ``classify_ode`` decided
    natively: a dict of the matching hints (possibly empty), or ``_MISSING``."""
    from symox.core.sympify import sympify
    from symox.polys.polytools import cancel
    native = getattr(_core, 'ode_power_series2', None)
    func = f(x)
    if native is None or _second_args(reduced_eq, func) is None or not isinstance(x, _core.Expr):
        return _MISSING
    try:
        at = sympify(point)
        if not isinstance(at, _core.Expr):
            return _MISSING
        value = native(reduced_eq, func, x, at, cancel)
    except Exception:
        return _MISSING
    if value is None:
        return _MISSING
    a3, b3, c3 = wilds
    if value[0] == 'ordinary':
        r = {a3: value[1], b3: value[2], c3: value[3]}
        r.update({'a3': a3, 'b3': b3, 'c3': c3, 'x0': point, 'terms': terms})
        return {"2nd_power_series_ordinary": r}
    if value[0] == 'regular':
        return {"2nd_power_series_regular": {'p': value[1], 'q': value[2], 'x0': point, 'terms': terms}}
    return {}


def series2_hook(matching_hints, reduced_eq, f, x, point, wilds, terms):
    """``classify_ode`` hook: ``True`` when the order-two power series hints
    were decided natively (and added to ``matching_hints``)."""
    native = power_series2(reduced_eq, f, x, point, wilds, terms)
    if native is _MISSING:
        return False
    matching_hints.update(native)
    return True


def nth_order_reducible(solver):
    """``NthOrderReducible._matches`` decided natively: ``(True, smallest)``,
    ``(False, None)``, or ``None`` when the reference matcher must run."""
    native = getattr(_core, 'ode_nth_order_reducible', None)
    problem = solver.ode_problem
    eq = problem.eq_preprocessed
    x = _linear_args(eq, problem.func)
    if native is None or x is None or x != problem.sym:
        return None
    try:
        value = native(eq, problem.func, x)
    except Exception:
        return None
    if value is None:
        return None
    if value is False:
        return False, None
    return True, value


def _expr_args(*values):
    return all(isinstance(v, _core.Expr) for v in values)


def nth_algebraic(eq, func, var, diffx, solve, simplify):
    """``NthAlgebraic._matches`` decided natively: the list of solution
    equations (empty for no match), or ``None`` when the reference body must
    run. ``diffx`` is the invertible derivative class of ``var``."""
    native = getattr(_core, 'ode_nth_algebraic', None)
    if native is None or not _expr_args(eq, func, var):
        return None
    template = diffx(var)
    if not isinstance(template, _core.Expr):
        return None
    try:
        value = native(eq, func, var, template, lambda e, s: solve(e, s, simplify=False), simplify)
    except Exception:
        return None
    return value


def factorable_match(eq, func, x, cancel, factor, solve):
    """``Factorable._matches`` decided natively: ``(matches, eqs)``, or
    ``None`` when the reference body must run."""
    native = getattr(_core, 'ode_factorable_match', None)
    if native is None or not _expr_args(eq, func, x):
        return None
    try:
        value = native(eq, func, x, cancel, factor, solve)
    except Exception:
        return None
    if value is None:
        return None
    if value is False:
        return False, None
    return True, value


def factorable_solutions(eqs, func, dsolve):
    """``Factorable._get_general_solution`` computed natively: ``(sols,
    unsolved)`` where ``unsolved`` is the last equation when no factor could
    be solved, or ``None`` when the reference body must run. An error raised
    by ``dsolve`` propagates, as in the reference."""
    native = getattr(_core, 'ode_factorable_solutions', None)
    if native is None or not eqs or not _expr_args(func, *eqs):
        return None
    return native(list(eqs), func, dsolve)


def bessel_match(solver, eq, factor, simplify, cancel, expand):
    """``SecondLinearBessel._matches`` decided natively: ``(True, rn)``,
    ``(False, None)``, or ``None`` when the reference matcher must run.
    The wild matches call the public ``Expr.match``."""
    native = getattr(_core, 'ode_bessel_match', None)
    problem = solver.ode_problem
    order = problem.order
    x = _second_args(eq, problem.func)
    if native is None or x is None or x != problem.sym or not isinstance(order, int) or order < 0:
        return None
    try:
        value = native(eq, problem.func, x, order, factor, simplify, cancel, expand, lambda e, p: e.match(p))
    except Exception:
        return None
    if value is None:
        return None
    if value is False:
        return False, None
    n, a4, c4, d4, b4 = value
    return True, {'n': n, 'a4': a4, 'd4': d4, 'c4': c4, 'b4': b4}


def bessel_solution(rn, func, x, C1, C2):
    """``SecondLinearBessel._get_general_solution`` computed natively: the
    equation, or ``None`` when the reference body must run."""
    native = getattr(_core, 'ode_bessel_solution', None)
    if native is None or not isinstance(rn, dict) or not _expr_args(func, x, C1, C2):
        return None
    try:
        params = tuple(rn[k] for k in ('n', 'a4', 'c4', 'd4', 'b4'))
    except KeyError:
        return None
    if not _expr_args(*params):
        return None
    try:
        return native(params, func, x, C1, C2)
    except Exception:
        return None
