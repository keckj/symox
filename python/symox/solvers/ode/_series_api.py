"""Exact local recurrences and first-order system inspection.

Series routines are bounded to 64 coefficients. Resonant Frobenius cases
requiring logarithmic reduction are reported explicitly.
"""
import re


def _S():
    import symox
    return symox


def _expressions(eq):
    S = _S()
    if isinstance(eq, (list, tuple, set)):
        return [_expressions(item) for item in eq]
    return S.sympify(eq)


def iter_numbered_constants(eq, start=1, prefix='C'):
    from .._differential import nodes
    S = _S()
    values = list(eq) if isinstance(eq, (list, tuple, set)) else [eq]
    names = set()
    for value in values:
        value = S.sympify(value)
        names.update(symbol.name for symbol in value.free_symbols)
        names.update(str(node.func) for node in nodes(value)
                     if getattr(node, 'is_UndefinedFunction', False))
    index = int(start)
    while True:
        name = prefix + str(index)
        index += 1
        if name not in names:
            yield S.Symbol(name)


def get_numbered_constants(eq, num=1, start=1, prefix='C'):
    if int(num) < 1:
        raise ValueError('num must be positive')
    stream = iter_numbered_constants(eq, start, prefix)
    result = tuple(next(stream) for _ in range(int(num)))
    return result[0] if num == 1 else result


def constant_renumber(expr, variables=None, newconstants=None):
    S = _S()
    values = list(expr) if isinstance(expr, (list, tuple, set)) else [expr]
    symbols = set().union(*(S.sympify(value).free_symbols for value in values))
    constants = (symbols-set(variables) if variables is not None else
                 {s for s in symbols if re.fullmatch(r'C\d+', s.name)})
    constants = sorted(constants, key=S.default_sort_key)
    protected = symbols-set(constants)
    stream = iter(newconstants) if newconstants is not None else iter_numbered_constants(list(protected))
    replacements = {constant: next(stream) for constant in constants}
    result = [S.sympify(value).xreplace(replacements) for value in values]
    if isinstance(expr, set):
        return set(result)
    if isinstance(expr, tuple):
        return tuple(result)
    return result if isinstance(expr, list) else result[0]


def solve_ics(sols, funcs, constants, ics):
    S = _S()
    from .._differential import nodes
    explicit = {}
    for solution in sols:
        for func in funcs:
            if solution.lhs == func:
                explicit[func] = solution.rhs
            elif solution.rhs == func:
                explicit[func] = solution.lhs
    equations = []
    for condition, value in ics.items():
        condition = S.sympify(condition)
        matched = False
        for func, rhs in explicit.items():
            x = func.args[0]
            if (getattr(condition, 'is_UndefinedFunction', False)
                    and condition.func == func.func and len(condition.args) == 1):
                evaluated = rhs.removeO().subs(x, condition.args[0])
                matched = True
            elif getattr(condition, 'is_Subs', False) or type(condition).__name__ == 'Subs':
                inner = condition.expr
                if inner.is_Derivative and inner.expr == func:
                    evaluated = rhs.removeO()
                    for variable, count in inner.variable_count:
                        evaluated = S.diff(evaluated, variable, int(count))
                    evaluated = evaluated.subs(dict(zip(condition.variables, condition.point)))
                    matched = True
            if matched:
                equations.append(evaluated-value)
                break
        if not matched:
            raise NotImplementedError('Initial condition must evaluate a solved function or its derivative')
    solutions = S.solve(equations, list(constants), dict=True)
    if len(solutions) != 1:
        raise ValueError('Initial conditions do not determine a unique set of constants')
    return solutions[0]


def ode_sol_simplicity(sol, func, trysolving=True):
    S = _S()
    if isinstance(sol, (tuple, list)):
        scores = [ode_sol_simplicity(item, func, trysolving) for item in sol]
        return S.oo if S.oo in scores else max(scores, default=-2)
    if sol.has(S.Integral):
        return S.oo
    if sol.lhs == func and not sol.rhs.has(func) or sol.rhs == func and not sol.lhs.has(func):
        return -2
    if trysolving:
        try:
            if S.solve(sol, func):
                return -1
        except (ValueError, NotImplementedError):
            pass
    return len(str(sol))


def odesimp(ode, eq, func, hint):
    S = _S()
    if isinstance(eq, (list, tuple)):
        return [odesimp(ode, item, func, hint) for item in eq]
    if eq.lhs == func:
        return S.Eq(func, S.simplify(eq.rhs))
    roots = S.solve(eq, func)
    result = [S.Eq(func, S.simplify(root)) for root in roots]
    return result[0] if len(result) == 1 else result or eq


def _parameters(match):
    S = _S()
    count = int(match.get('terms', 6))
    if not 2 <= count <= 64:
        raise ValueError('Series routines require 2 to 64 terms')
    return S.sympify(match.get('x0', match.get('f0', 0))), count


def _finite(value):
    S = _S()
    return not value.has(S.oo, -S.oo, S.zoo, S.nan)


def _slope(eq, func):
    from ._native_ode import _data
    S = _S()
    expression, func, x, order, coordinates = _data(eq, func)
    if order != 1:
        raise ValueError('Expected a first-order differential equation')
    roots = S.solve(expression, coordinates[1])
    if len(roots) != 1:
        raise NotImplementedError('The highest derivative must have a unique explicit branch')
    return roots[0]


def ode_1st_power_series(eq, func, order, match):
    S = _S()
    if order != 1:
        raise ValueError('Expected order one')
    x, y = func.args[0], S.Dummy('dependent_value')
    point, count = _parameters(match)
    initial = S.sympify(match.get('f0val', get_numbered_constants(eq)))
    slope = _slope(eq, func).xreplace({func: y})
    derivative, result = slope, initial
    for k in range(1, count):
        value = S.simplify(derivative.subs({x: point, y: initial}))
        if not _finite(value):
            raise ValueError('The slope is not analytic at the initial point')
        result += value*(x-point)**k/S.factorial(k)
        derivative = S.diff(derivative, x)+slope*S.diff(derivative, y)
    return S.Eq(func, result+S.Order((x-point)**count, (x, point)))


def _linear_second(eq, func):
    from ._native_ode import _data
    from .._differential import linear_data
    S = _S()
    expression, func, x, order, coordinates = _data(eq, func)
    data = linear_data(expression, coordinates)
    if order != 2 or data is None:
        raise ValueError('Expected a linear second-order differential equation')
    coefficients, rest = data
    return x, [S.cancel(value/coefficients[2]) for value in (coefficients[1], coefficients[0], rest)]


def _taylor(value, x, point, h, count):
    S = _S()
    polynomial = S.series(value.subs(x, h+point), h, 0, count).removeO().expand()
    if not _finite(polynomial) or any(polynomial.coeff(h, k).has(h) for k in range(count)):
        raise ValueError('Nonanalytic coefficient at expansion point')
    if S.simplify(polynomial-sum((polynomial.coeff(h, k)*h**k for k in range(count)), S.S.Zero)) != 0:
        raise ValueError('Expected an ordinary Taylor expansion')
    return polynomial


def ode_2nd_power_series_ordinary(eq, func, order, match):
    S = _S()
    x, (p, q, r) = _linear_second(eq, func)
    point, count = _parameters(match)
    h = S.Dummy('local_coordinate')
    p, q, r = [_taylor(value, x, point, h, count) for value in (p, q, r)]
    c0, c1 = get_numbered_constants(eq, 2)
    polynomial = c0+c1*h
    for k in range(count-2):
        coefficient = (-p*polynomial.diff(h)-q*polynomial-r).expand().coeff(h, k)
        polynomial += S.cancel(coefficient/((k+2)*(k+1)))*h**(k+2)
    return S.Eq(func, polynomial.subs(h, x-point)+S.Order((x-point)**count, (x, point)))


def ode_2nd_power_series_regular(eq, func, order, match):
    S = _S()
    x, (p, q, forcing) = _linear_second(eq, func)
    if forcing != 0:
        raise NotImplementedError('Frobenius expansion requires a homogeneous equation')
    point, count = _parameters(match)
    h, r = S.Dummy('local_coordinate'), S.Dummy('indicial_root')
    p = _taylor((x-point)*p, x, point, h, count)
    q = _taylor((x-point)**2*q, x, point, h, count)
    p0, q0 = p.coeff(h, 0), q.coeff(h, 0)
    indicial = r*(r-1)+p0*r+q0
    roots = S.solve(indicial, r)
    if len(roots) != 2:
        raise NotImplementedError('Repeated indicial roots require logarithmic reduction')
    solutions = []
    for root in roots:
        coefficients = [S.S.One]
        for n in range(1, count):
            numerator = sum((((root+k)*p.coeff(h, n-k)+q.coeff(h, n-k))*coefficients[k]
                             for k in range(n)), S.S.Zero)
            denominator = S.simplify(indicial.subs(r, root+n))
            if denominator == 0:
                if S.simplify(numerator) != 0:
                    raise NotImplementedError('Resonant indicial roots require a logarithmic solution')
                coefficients.append(S.S.Zero)
            else:
                coefficients.append(S.cancel(-numerator/denominator))
        solutions.append((x-point)**root*(sum((a*(x-point)**k for k, a in enumerate(coefficients)), S.S.Zero)
                                         +S.Order((x-point)**count, (x, point))))
    c1, c2 = get_numbered_constants(eq, 2)
    return S.Eq(func, c1*solutions[0]+c2*solutions[1])


def classify_sysode(eq, funcs=None, **kwargs):
    from .._differential import equation, nodes, linear_data
    from ..deutils import ode_order
    S = _S()
    if kwargs:
        raise NotImplementedError('System classification options are not implemented')
    equations = [equation(item) for item in eq]
    if funcs is None:
        funcs = sorted({node.expr for expr in equations for node in nodes(expr)
                        if node.is_Derivative and getattr(node.expr, 'is_UndefinedFunction', False)}, key=S.default_sort_key)
    funcs = list(funcs)
    if not funcs or len(equations) != len(funcs):
        raise ValueError('Supply one equation per dependent function')
    t = funcs[0].args[0]
    orders = {func: max(ode_order(expr, func) for expr in equations) for func in funcs}
    coordinates = [S.diff(func, t, k) for func in funcs for k in range(orders[func]+1)]
    linear = all(linear_data(expr, coordinates) is not None for expr in equations)
    coefficients = {}
    placeholders = [S.Dummy('coordinate') for _ in coordinates]
    for i, expr in enumerate(equations):
        expanded = expr.xreplace(dict(zip(coordinates, placeholders))).expand()
        for func in funcs:
            for k in range(orders[func]+1):
                coordinate = S.diff(func, t, k)
                coefficients[i, func, k] = expanded.coeff(placeholders[coordinates.index(coordinate)]).xreplace(dict(zip(placeholders, coordinates)))
    match = {'eq': equations, 'func': funcs, 'order': orders, 'no_of_equation': len(equations),
             'is_linear': linear, 'func_coeff': coefficients, 'type_of_equation': None}
    if len(set(orders.values())) == 1:
        order = next(iter(orders.values()))
        checker = {(True, 2, 1): check_linear_2eq_order1,
                   (False, 2, 1): check_nonlinear_2eq_order1,
                   (False, 3, 1): check_nonlinear_3eq_order1,
                   (False, 2, 2): check_nonlinear_2eq_order2,
                   (False, 3, 2): check_nonlinear_3eq_order2}.get((linear, len(funcs), order))
        if checker:
            match['type_of_equation'] = checker(equations, funcs, coefficients)
    return match


def _slopes(eq, funcs):
    S = _S()
    from ._native_systems import canonical_odes
    branches = canonical_odes(eq, funcs, funcs[0].args[0])
    if len(branches) != 1:
        raise NotImplementedError('System must have a unique derivative branch')
    explicit = branches[0]
    return [next(item.rhs for item in explicit if item.lhs == S.diff(func, func.args[0])) for func in funcs]


def check_linear_2eq_order1(eq, func, func_coef):
    S = _S()
    t = func[0].args[0]
    try:
        slopes = _slopes(eq, func)
    except (ValueError, NotImplementedError):
        return None
    matrix = S.Matrix([[slope.diff(f) for f in func] for slope in slopes])
    if any(S.simplify(slope-sum(matrix[i, j]*func[j] for j in range(2))) != 0 for i, slope in enumerate(slopes)):
        return None
    if all(not entry.has(t) for entry in matrix):
        return None
    a = S.Dummy('constant_ratio')
    roots = S.solve(matrix[1, 0]+(matrix[1, 1]-matrix[0, 0])*a-matrix[0, 1]*a*a, a)
    return 'type6' if any(not root.has(t) for root in roots) else 'type7'


def check_nonlinear_2eq_order1(eq, func, func_coef):
    S = _S()
    t = func[0].args[0]
    try:
        slopes = _slopes(eq, func)
    except (ValueError, NotImplementedError):
        return None
    u, v = S.Dummy('u'), S.Dummy('v')
    converted = [s.xreplace(dict(zip(func, (u, v)))) for s in slopes]
    ratio = S.cancel(converted[0]/converted[1])
    separated = S.separatevars(ratio, symbols=[u, v], dict=True)
    autonomous = all(not slope.has(t) for slope in converted)
    if separated is not None and not ratio.has(t):
        if autonomous:
            upper = separated[u]
            if upper == u or (upper.is_Pow and upper.base == u) or not upper.has(u):
                return 'type1'
            if getattr(upper.func, '__name__', '') == 'exp':
                return 'type2'
        return 'type4'
    return 'type3' if autonomous else None


def check_nonlinear_3eq_order1(eq, func, func_coef):
    S = _S()
    try:
        slopes = _slopes(eq, func)
    except (ValueError, NotImplementedError):
        return None
    t = func[0].args[0]
    factors = [S.cancel(slopes[i]/(func[(i+1)%3]*func[(i+2)%3])) for i in range(3)]
    if any(value == 0 for value in factors):
        return None
    ratios = [S.cancel(value/factors[0]) for value in factors]
    if any(value.has(t, *func) for value in ratios):
        return None
    if all(not value.has(t, *func) for value in factors):
        a, b, c = factors
        return 'type1' if S.simplify(a+b+c+a*b*c) == 0 else None
    return 'type2'


def check_nonlinear_2eq_order2(eq, func, func_coef):
    """No second-order nonlinear system family is classified (as upstream)."""
    return None


def check_nonlinear_3eq_order2(eq, func, func_coef):
    """No second-order nonlinear system family is classified (as upstream)."""
    return None


def sysode_linear_2eq_order1(match_):
    from ._native_systems import dsolve_system
    return dsolve_system(match_['eq'], match_['func'])


def _nonlinear_system(match_):
    """Solve triangular first-order systems, substituting each exact solution."""
    S = _S()
    from .._differential import replace_function
    funcs, equations = list(match_['func']), list(match_['eq'])
    results = []
    while funcs:
        for i, expr in enumerate(equations):
            candidates = [func for func in funcs if expr.has(func)]
            if len(candidates) != 1:
                continue
            func = candidates[0]
            solution = S.dsolve(expr, func)
            if isinstance(solution, (list, tuple)) or solution.lhs != func or solution.rhs.has(func):
                raise NotImplementedError('Triangular system requires a unique explicit scalar solution')
            occupied = set().union(*(item.free_symbols for item in results)) if results else set()
            constants = {s for s in solution.free_symbols-expr.free_symbols if re.fullmatch(r'C\d+', s.name)}
            stream = iter_numbered_constants(list(occupied)+[expr])
            solution = solution.xreplace({constant: next(stream) for constant in sorted(constants, key=S.default_sort_key)})
            results.append(solution)
            funcs.remove(func)
            equations.pop(i)
            equations = [S.simplify(replace_function(item, func, solution.rhs)) for item in equations]
            break
        else:
            raise NotImplementedError('Nonlinear system is coupled; triangular reduction is unavailable')
    return results


def sysode_nonlinear_2eq_order1(match_):
    return _nonlinear_system(match_)


def sysode_nonlinear_3eq_order1(match_):
    return _nonlinear_system(match_)
