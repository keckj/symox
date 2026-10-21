"""Recognize implemented scalar ODE families using exact symbolic tests."""
from .._differential import equation, dependent, nodes, linear_data, fresh_constant
from ..deutils import ode_order
from ._constants import constantsimp
from .lie_group import infinitesimals, checkinfsol, solve_lie

allhints = ('separable', '1st_exact', '1st_linear', 'Bernoulli',
            'lie_group',
            'nth_linear_constant_coeff_homogeneous',
            'nth_linear_constant_coeff_undetermined_coefficients',
            'nth_linear_constant_coeff_variation_of_parameters',
            'nth_linear_euler_eq_homogeneous', '1st_power_series',
            '2nd_power_series_ordinary', '2nd_power_series_regular')


def _zero(value):
    from ... import simplify, cancel
    value = simplify(cancel(value))
    return value == 0 or value.is_zero is True


def _data(eq, func=None):
    from ... import Derivative
    expr = equation(eq)
    func = dependent(expr, func)
    order = ode_order(expr, func)
    if order > 256:
        raise NotImplementedError('Differential order exceeds 256')
    x = func.args[0]
    coordinates = [func]+[Derivative(func, (x, k)) for k in range(1, order+1)]
    if any(node.is_Derivative and node.has(func) and node not in coordinates for node in nodes(expr)):
        raise NotImplementedError('Classification requires derivatives of the dependent function itself')
    return expr, func, x, order, coordinates


def classify_ode(eq, func=None, dict=False, ics=None, *, prep=True, xi=None, eta=None, n=None, **kwargs):
    from ... import Dummy, diff, cancel, expand, separatevars, Add
    if kwargs:
        raise NotImplementedError('These ODE classification options are not implemented')
    if (xi is None)!=(eta is None):
        raise ValueError('Supply both xi and eta for a custom Lie generator')
    expr, func, x, order, coordinates = _data(eq, func)
    matches = {}
    linear = linear_data(expr, coordinates) if order else None
    if linear is not None:
        coefficients, rest = linear
        # Other occurrences of the dependent function, including derivatives
        # inside coefficient functions, are never independent coefficients.
        if any(value.has(func) for value in coefficients+[rest]):
            linear = None
        else:
            normalized = [cancel(value/coefficients[-1]) for value in coefficients]
            match = {'coefficients': tuple(coefficients), 'forcing': -rest}
            if order == 1:
                matches['1st_linear'] = match
                matches['Bernoulli'] = match
            if all(not value.has(x) for value in normalized):
                if _zero(rest):
                    matches['nth_linear_constant_coeff_homogeneous'] = match
                else:
                    # Variation of parameters is supported by the native
                    # second-order solver. Undetermined coefficients only
                    # accepts a finite polynomial/exponential/trig basis.
                    if order == 2:
                        matches['nth_linear_constant_coeff_variation_of_parameters'] = match
                    if _finite_forcing(cancel(-rest/coefficients[-1]), x):
                        matches['nth_linear_constant_coeff_undetermined_coefficients'] = match
            if _zero(rest) and all(not cancel(value*x**(order-k)/coefficients[-1]).has(x)
                                   for k, value in enumerate(coefficients)):
                matches['nth_linear_euler_eq_homogeneous'] = match
    if order == 1:
        y, yp = Dummy('dependent_value'), Dummy('dependent_derivative')
        polynomial = expand(cancel(expr.xreplace({func:y, coordinates[1]:yp})).as_numer_denom()[0])
        a = polynomial.coeff(yp)
        b = expand(polynomial-a*yp)
        if not a.has(yp) and not b.has(yp) and not _zero(a):
            matches['lie_group'] = {'xi':xi,'eta':eta}
            if xi is not None:
                from ... import Function
                generator = {Function('xi')(x,func):xi,Function('eta')(x,func):eta}
                if not checkinfsol(expr,[generator],func)[0][0]:
                    raise ValueError('The supplied infinitesimals are not a symmetry of this ODE')
            slope = cancel(-b/a)
            separated = separatevars(slope, symbols=[x, y], dict=True)
            if separated is not None:
                matches['separable'] = {'x':separated[x], 'y':separated[y],
                                        'coefficient':separated['coeff']}
            if _zero(diff(b, y)-diff(a, x)):
                matches['1st_exact'] = {'M':b, 'N':a}
            if 'Bernoulli' not in matches:
                parts = []
                for term in Add.make_args(expand(slope)):
                    separated_term = separatevars(term, symbols=[x, y], dict=True)
                    if separated_term is None:
                        break
                    factor = separated_term[y]
                    if factor == y:
                        power = 1
                    elif factor.is_Pow and factor.base == y and not factor.exp.has(x, y):
                        power = factor.exp
                    else:
                        break
                    parts.append(power)
                else:
                    if parts and len(set(parts)-{1}) == 1:
                        matches['Bernoulli'] = {'slope':slope}
    if n is not None:
        from ._series_api import _linear_second, _taylor
        from ... import S, Dummy
        count = int(n)
        if not 2 <= count <= 64:raise ValueError('Series classification requires 2 to 64 terms')
        if order == 1 and 'lie_group' in matches:
            matches['1st_power_series'] = {'terms': count, 'f0': S.Zero}
        if order == 2 and linear is not None:
            h = Dummy('local_coordinate')
            _, (p, q, rest) = _linear_second(expr,func)
            try:
                for value in (p,q,rest):_taylor(value,x,0,h,count)
                matches['2nd_power_series_ordinary'] = {'terms':count,'x0':S.Zero}
            except (ValueError,NotImplementedError,ZeroDivisionError):
                if _zero(rest):
                    try:
                        _taylor(x*p,x,0,h,count);_taylor(x*x*q,x,0,h,count)
                        matches['2nd_power_series_regular'] = {'terms':count,'x0':S.Zero}
                    except (ValueError,NotImplementedError,ZeroDivisionError):pass
    hints = tuple(hint for hint in allhints if hint in matches)
    if not dict:
        return hints
    return {'order':order, **matches, 'default':hints[0] if hints else None, 'ordered_hints':hints}


def _finite_forcing(value, x):
    from ... import Add, Mul, diff
    for term in Add.make_args(value.expand()):
        for factor in Mul.make_args(term):
            if not factor.has(x) or factor == x:
                continue
            if factor.is_Pow and factor.base == x and factor.exp.is_integer is True and factor.exp.is_nonnegative is True:
                if int(factor.exp) <= 256:
                    continue
            if getattr(factor.func, '__name__', '') in ('exp', 'sin', 'cos'):
                if not diff(factor.args[0], x).has(x):
                    continue
            return False
    return True


def dsolve(eq, func=None, hint='default', simplify=True, ics=None, **kwargs):
    from ... import _core, Eq, Dummy, Derivative, integrate, diff, solve, cancel
    if hint in ('1st_power_series', '2nd_power_series_ordinary', '2nd_power_series_regular'):
        point, count = kwargs.pop('x0', 0), kwargs.pop('n', 6)
        if kwargs:
            raise NotImplementedError('Unsupported series options')
        expr, selected, x, order, coordinates = _data(eq, func)
        match = {'x0': point, 'terms': count}
        if ics and hint == '1st_power_series':
            key = selected.func(point)
            if set(ics) != {key}:
                raise NotImplementedError('First-order series require one value at x0')
            match['f0val'] = ics[key]
        solver = {'1st_power_series': ode_1st_power_series,
                  '2nd_power_series_ordinary': ode_2nd_power_series_ordinary,
                  '2nd_power_series_regular': ode_2nd_power_series_regular}[hint]
        result = solver(expr, selected, order, match)
        if ics and order == 2:
            constants = result.free_symbols-expr.free_symbols-{x}
            result = result.subs(solve_ics([result], [selected], constants, ics))
        return result
    xi,eta = kwargs.pop('xi',None),kwargs.pop('eta',None)
    if not simplify or kwargs:
        raise NotImplementedError('These dsolve options are not implemented')
    if (xi is None)!=(eta is None):
        raise ValueError('Supply both xi and eta for a custom Lie generator')
    if isinstance(eq,(list,tuple)):
        if hint != 'default':
            raise NotImplementedError('Named hints for ODE systems are not implemented')
        from ._native_systems import dsolve_system, ODENonlinearError
        try:return dsolve_system(eq,func,ics=ics)
        except ODENonlinearError:
            match = classify_sysode(eq,func)
            if len(match['func']) not in (2,3):raise
            result = (sysode_nonlinear_2eq_order1(match) if len(match['func']) == 2 else
                      sysode_nonlinear_3eq_order1(match))
            if ics:
                symbols = set().union(*(item.free_symbols for item in result))
                protected = set().union(*(item.free_symbols for item in match['eq']))
                substitutions = solve_ics(result, match['func'], symbols-protected, ics)
                result = [item.subs(substitutions) for item in result]
            return result
    if hint == 'default':
        try:
            return _core.dsolve(eq, func, ics=ics)
        except NotImplementedError as original:
            try:
                expr, selected, x, order, coordinates = _data(eq, func)
                hints = classify_ode(expr, selected)
            except (ValueError, TypeError, NotImplementedError):
                raise original
            if '1st_exact' in hints:
                return dsolve(expr, selected, hint='1st_exact', ics=ics)
            linear = linear_data(expr, coordinates)
            if linear is not None and not _zero(linear[0][-1]-1):
                try:
                    return _core.dsolve(cancel(expr/linear[0][-1]), selected, ics=ics)
                except NotImplementedError:
                    pass
            if 'lie_group' in hints:
                return solve_lie(expr,selected,xi=xi,eta=eta,ics=ics)
            raise original
    if hint not in allhints:
        raise NotImplementedError('Unknown or unimplemented ODE hint: '+str(hint))
    if hint == 'lie_group':
        return solve_lie(eq,func,xi=xi,eta=eta,ics=ics)
    expr, func, x, order, coordinates = _data(eq, func)
    if hint not in classify_ode(expr, func):
        raise ValueError('The equation does not match hint '+hint)
    if hint != '1st_exact':
        return _core.dsolve(expr, func, ics=ics)
    if ics:
        raise NotImplementedError('Initial conditions for implicit exact equations are not implemented')
    y, yp = Dummy('dependent_value'), Dummy('dependent_derivative')
    reduced = cancel(expr.xreplace({func:y, coordinates[1]:yp})).as_numer_denom()[0].expand()
    a = reduced.coeff(yp)
    b = (reduced-a*yp).expand()
    potential = integrate(b, x)
    potential += integrate(cancel(a-diff(potential, y)), y)
    constant = fresh_constant(expr)
    try:
        roots = solve(potential-constant, y)
    except NotImplementedError:
        roots = []
    if roots:
        solutions = [Eq(func, root) for root in roots]
        return solutions[0] if len(solutions) == 1 else solutions
    return Eq(potential.xreplace({y:func}), constant)


def checkodesol(ode, sol, func=None, order='auto', solve_for_func=True):
    from ... import _core, Equality, Eq, Dummy, diff, cancel, simplify
    expression = equation(ode)
    if isinstance(sol, (list, tuple, set)):
        return [checkodesol(expression, item, func, order, solve_for_func) for item in sol]
    if not isinstance(sol, Equality):
        return _core.checkodesol(expression, sol, func)
    # Preserve the native explicit checker, including its conservative
    # behavior on held derivatives of composite expressions.
    if func is None:
        try:
            selected = dependent(expression)
        except ValueError:
            return _core.checkodesol(expression, sol, func)
    else:
        selected = func
    if sol.lhs == selected:
        return _core.checkodesol(expression, sol, selected)
    if sol.rhs == selected:
        return _core.checkodesol(expression, Eq(selected, sol.lhs), selected)
    if (not getattr(selected, 'is_UndefinedFunction', False) or len(selected.args) != 1
            or not selected.args[0].is_Symbol):
        return _core.checkodesol(expression, sol, selected)
    expr, selected, x, actual_order, coordinates = _data(expression, selected)
    if order != 'auto' and order != actual_order:
        raise ValueError('The supplied ODE order does not match the equation')
    y = Dummy('implicit_dependent_coordinate')
    relation = (sol.lhs-sol.rhs).xreplace({selected:y})
    denominator = diff(relation, y)
    if _zero(denominator):
        raise ValueError('The proposed solution does not constrain the dependent function')
    first = cancel(-diff(relation, x)/denominator)
    replacements = {selected:y}
    derivative = first
    for index, coordinate in enumerate(coordinates[1:]):
        if index:
            derivative = cancel(diff(derivative, x)+diff(derivative, y)*first)
        replacements[coordinate] = derivative
    residual = simplify(cancel(expression.xreplace(replacements)))
    return _zero(residual), residual.xreplace({y:selected})


__all__ = ['classify_ode', 'dsolve', 'checkodesol', 'allhints']

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox', 'homogeneous_order')

from ._series_api import (get_numbered_constants, iter_numbered_constants, constant_renumber,
    solve_ics, ode_sol_simplicity, odesimp, ode_1st_power_series,
    ode_2nd_power_series_ordinary, ode_2nd_power_series_regular, classify_sysode,
    check_linear_2eq_order1, check_nonlinear_2eq_order1, check_nonlinear_2eq_order2,
    check_nonlinear_3eq_order1, check_nonlinear_3eq_order2, sysode_linear_2eq_order1,
    sysode_nonlinear_2eq_order1, sysode_nonlinear_3eq_order1)
__all__ += ['get_numbered_constants', 'iter_numbered_constants', 'constant_renumber',
    'solve_ics', 'ode_sol_simplicity', 'odesimp', 'ode_1st_power_series',
    'ode_2nd_power_series_ordinary', 'ode_2nd_power_series_regular', 'classify_sysode',
    'check_linear_2eq_order1', 'check_nonlinear_2eq_order1', 'check_nonlinear_2eq_order2',
    'check_nonlinear_3eq_order1', 'check_nonlinear_3eq_order2', 'sysode_linear_2eq_order1',
    'sysode_nonlinear_2eq_order1', 'sysode_nonlinear_3eq_order1']
