"""First-order linear PDEs and variable separation on native expressions."""
from ._differential import equation, dependent, linear_data, replace_function, nodes
from .deutils import ode_order

allhints = ('1st_linear_constant_coeff_homogeneous', '1st_linear_constant_coeff',
            '1st_linear_constant_coeff_Integral', '1st_linear_variable_coeff')


def _zero(value):
    from .. import simplify, cancel
    value = simplify(cancel(value))
    return value == 0 or value.is_zero is True


def _data(eq, func=None):
    from .. import Derivative
    expression = equation(eq)
    func = dependent(expression, func, dimension=2)
    order = ode_order(expression, func)
    if order > 256:
        raise NotImplementedError('Differential order exceeds 256')
    x, y = func.args
    coordinates = [Derivative(func, x), Derivative(func, y), func]
    linear = linear_data(expression, coordinates) if order == 1 else None
    if linear is not None:
        coefficients, rest = linear
        if any(value.has(func) for value in coefficients+[rest]):
            linear = None
        elif _zero(coefficients[0]) and _zero(coefficients[1]):
            linear = None
    return expression, func, order, linear


def classify_pde(eq, func=None, dict=False, *, prep=True, **kwargs):
    if kwargs:
        raise NotImplementedError('These PDE classification options are not implemented')
    expression, func, order, linear = _data(eq, func)
    matches = {}
    if linear is not None:
        coefficients, rest = linear
        match = {'coefficients':tuple(coefficients), 'forcing':-rest}
        if not any(value.has(*func.args) for value in coefficients):
            if _zero(rest):
                matches[allhints[0]] = match
            else:
                matches[allhints[1]] = match
                matches[allhints[2]] = match
        else:
            matches[allhints[3]] = match
    hints = tuple(hint for hint in allhints if hint in matches)
    return ({'order':order, **matches, 'default':hints[0] if hints else None,
             'ordered_hints':hints} if dict else hints)


def pdsolve(eq, func=None, hint='default', dict=False, solvefun=None, **kwargs):
    from .. import Eq, Dummy, Function, exp, integrate, Integral, Subs, cancel, simplify
    if kwargs:
        raise NotImplementedError('These pdsolve options are not implemented')
    expression, func, order, linear = _data(eq, func)
    hints = classify_pde(expression, func)
    if not hints:
        raise NotImplementedError('No implemented first-order linear PDE method matches')
    if hint in ('all', 'all_Integral'):
        chosen = [item for item in hints if hint != 'all_Integral' or item.endswith('_Integral')
                  or item == allhints[0] or item == allhints[3]]
        result = {'order':order, 'default':hints[0]}
        for item in chosen:
            try:
                result[item] = pdsolve(expression, func, hint=item, solvefun=solvefun)
            except NotImplementedError as error:
                result[item] = error
        return result
    if hint == 'default':
        hint = hints[0]
    elif hint not in allhints:
        raise ValueError('Unrecognized PDE hint: '+str(hint))
    elif hint not in hints:
        raise ValueError('The equation does not match PDE hint '+hint)
    if solvefun is None:
        occupied = {getattr(node.func, '__name__', '') for node in nodes(expression)
                    if node.is_UndefinedFunction}
        name = 'F'
        while name in occupied:
            name += '_'
        solvefun = Function(name)
    try:
        probe = solvefun(Dummy('characteristic_probe'))
    except TypeError as error:
        raise TypeError('solvefun must be an undefined function head') from error
    if not getattr(probe, 'is_UndefinedFunction', False):
        raise TypeError('solvefun must be an undefined function head')
    a, b, c = linear[0]
    source = -linear[1]
    x, y = func.args
    if hint == allhints[3]:
        value = _variable_solution(a, b, c, source, x, y, solvefun)
    else:
        norm = cancel(a*a+b*b)
        if _zero(norm):
            raise NotImplementedError('Isotropic complex characteristic directions are not implemented')
        t, invariant = Dummy('characteristic_time'), Dummy('characteristic_invariant')
        time_value = a*x+b*y
        invariant_value = b*x-a*y
        forcing = source.xreplace({x:(a*t+b*invariant)/norm, y:(b*t-a*invariant)/norm})
        mu = exp(c*t/norm)
        integrand = simplify(forcing*mu/norm)
        if _zero(source):
            particular = 0
        elif hint.endswith('_Integral'):
            # The one-bound integral is an antiderivative evaluated at the
            # characteristic coordinate; its dummy cannot capture user names.
            particular = Integral(integrand.xreplace({invariant:invariant_value}), (t, time_value))
        else:
            particular = integrate(integrand, t).xreplace({invariant:invariant_value})
            if particular.has(t):
                particular = Subs(particular, t, time_value).doit()
        value = exp(-c*time_value/norm)*(solvefun(invariant_value)+particular)
    result = Eq(func, value)
    return {func:value} if dict else result


def _variable_solution(a, b, c, source, x, y, arbitrary):
    from .. import integrate, exp, cancel, simplify, diff
    if _zero(a):
        mu = exp(integrate(cancel(c/b), y))
        return (arbitrary(x)+integrate(cancel(source*mu/b), y))/mu
    if _zero(b):
        mu = exp(integrate(cancel(c/a), x))
        return (arbitrary(y)+integrate(cancel(source*mu/a), x))/mu
    ax, by = cancel(a/x), cancel(b/y)
    if not ax.has(x, y) and not by.has(x, y) and not c.has(x, y) and _zero(source):
        # Diagonal dilation vector fields admit a simple algebraic invariant.
        return arbitrary(y*x**(-by/ax))*x**(-c/ax)
    if not a.has(y) and not b.has(x) and _zero(c) and _zero(source):
        return arbitrary(integrate(1/a, x)-integrate(1/b, y))
    # A shear with one constant velocity can be straightened explicitly.
    if not a.has(x, y) and not b.has(y) and _zero(c) and _zero(source):
        return arbitrary(y-integrate(cancel(b/a), x))
    if not b.has(x, y) and not a.has(x) and _zero(c) and _zero(source):
        return arbitrary(x-integrate(cancel(a/b), y))
    raise NotImplementedError('This variable-coefficient characteristic system is not implemented')


def checkpdesol(pde, sol, func=None, solve_for_func=True):
    from .. import Equality, simplify, cancel, solve
    expression = equation(pde)
    func = dependent(expression, func, dimension=2)
    if isinstance(sol, (list, tuple, set)):
        return [checkpdesol(expression, item, func, solve_for_func=solve_for_func) for item in sol]
    if isinstance(sol, Equality):
        if sol.lhs == func:
            value = sol.rhs
        elif sol.rhs == func:
            value = sol.lhs
        elif solve_for_func:
            values = solve(sol.lhs-sol.rhs, func)
            if len(values) != 1:
                raise NotImplementedError('Cannot isolate a unique PDE solution')
            value = values[0]
        else:
            raise NotImplementedError('Implicit PDE solution checks need function isolation')
    else:
        from ..calculus._analysis import scalar
        value = scalar(sol)
    residual = replace_function(expression, func, value)
    residual = simplify(cancel(_integral_linearity(residual)))
    return _zero(residual), residual


def _integral_linearity(expression):
    """Combine residual integrals with identical binders and endpoints.

    Coefficients may move inside only when they are independent of every
    bound variable. No endpoint substitutions or numerical integration are
    used to certify the resulting cancellation.
    """
    from .. import Add, Mul, Integral, simplify, S
    groups, remaining = {}, []
    for term in Add.make_args(expression.expand()):
        factors = Mul.make_args(term)
        integrals = [factor for factor in factors if factor.is_Integral]
        if len(integrals) != 1:
            remaining.append(term)
            continue
        integral = integrals[0]
        coefficient = Mul(*(factor for factor in factors if factor is not integral))
        limits = tuple(tuple(limit) for limit in integral.limits)
        if coefficient.has(*(limit[0] for limit in limits)):
            remaining.append(term)
            continue
        groups.setdefault(limits, []).append(coefficient*integral.function)
    for limits, terms in groups.items():
        integrand = simplify(Add(*terms))
        if not _zero(integrand):
            remaining.append(Integral(integrand, *limits))
    return Add(*remaining)


def pde_separate(eq, fun, sep, strategy='mul'):
    from .. import sympify, Add, Mul, expand, cancel, separatevars, S
    if strategy not in ('mul', 'add'):
        raise ValueError('Unknown separation strategy: '+str(strategy))
    expression = equation(eq)
    fun = sympify(fun)
    if not fun.is_UndefinedFunction or not all(arg.is_Symbol for arg in fun.args):
        raise ValueError('Separation requires an applied function of independent symbols')
    sep = tuple(sympify(value) for value in sep)
    if len(sep) < 2 or any(not value.is_UndefinedFunction for value in sep):
        raise ValueError('Supply at least two separated functions')
    arguments = [arg for value in sep for arg in value.args]
    if len(arguments) != len(fun.args):
        raise ValueError('Variable counts do not match')
    if len(set(arguments)) != len(arguments):
        raise ValueError('Duplicate substitution arguments detected')
    if set(arguments) != set(fun.args):
        raise ValueError('Arguments do not match')
    value = Add(*sep) if strategy == 'add' else Mul(*sep)
    expression = replace_function(expression, fun, value)
    if strategy == 'mul':
        expression = cancel(expression/value)
    expression = expand(expression)
    selected = set(sep[0].args)
    others = set(arguments)-selected
    # First divide out a common other-variable factor from terms carrying
    # derivatives in the selected group. This also separates nonlinear
    # exponential equations under an additive ansatz.
    divisors = []
    for term in Add.make_args(expression):
        if any(node.is_Derivative and not node.has(*others) and node.has(*selected)
               for node in nodes(term)):
            independent, dependent_part = term.as_independent(*selected, as_Add=False)
            if dependent_part.has(*others):
                return None
            divisors.append(independent)
    candidates = [S.One]+divisors
    for divisor in candidates:
        if _zero(divisor):
            continue
        normalized = expand(cancel(expression/divisor))
        left, mixed = S.Zero, []
        for term in Add.make_args(normalized):
            if not term.has(*others):
                left += term
            else:
                independent, part = term.as_independent(*selected, as_Add=False)
                if part.has(*others):
                    break
                mixed.append((independent, part))
        else:
            if not mixed or left == 0:
                continue
            common = mixed[0][1]
            if all(_zero(part-common) for _, part in mixed):
                lhs = cancel(left/common).expand()
                rhs = -sum((independent for independent, _ in mixed), S.Zero)
                if not lhs.has(*others) and not rhs.has(*selected):
                    return [lhs, rhs.expand()]
    return None


def pde_separate_add(eq, fun, sep):
    return pde_separate(eq, fun, sep, strategy='add')


def pde_separate_mul(eq, fun, sep):
    return pde_separate(eq, fun, sep, strategy='mul')


__all__ = ['pdsolve', 'classify_pde', 'checkpdesol', 'pde_separate',
           'pde_separate_add', 'pde_separate_mul']
