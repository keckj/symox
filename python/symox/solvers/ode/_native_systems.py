"""Exact linear ODE systems with constant or commuting coefficient matrices."""
from .._differential import equation, dependent, nodes, linear_data
from ..deutils import ode_order


class ODEOrderError(ValueError):
    pass


class ODENonlinearError(ValueError):
    pass


def _inputs(eqs, funcs, t):
    from ... import sympify
    eqs, funcs, t = [equation(eq) for eq in eqs], list(funcs), sympify(t)
    if not t.is_Symbol or not funcs or len(set(funcs)) != len(funcs):
        raise ValueError('Supply distinct dependent functions and an independent symbol')
    if len(eqs) > 64 or len(funcs) > 64:
        raise NotImplementedError('System inspection is limited to 64 equations and functions')
    funcs = [dependent(eqs[0], func) for func in funcs]
    if any(func.args != (t,) for func in funcs):
        raise ValueError('All system functions must depend on the same independent symbol')
    return eqs, funcs, t


def linear_ode_to_matrix(eqs, funcs, t, order):
    from ... import Matrix, Derivative, S
    eqs, funcs, t = _inputs(eqs, funcs, t)
    if not isinstance(order, int) or isinstance(order, bool) or order < 0:
        raise ValueError('System order must be a nonnegative integer')
    if order > 256:
        raise NotImplementedError('System order exceeds 256')
    if any(ode_order(eq, func) > order for eq in eqs for func in funcs):
        raise ODEOrderError('Cannot represent the system at the requested differential order')
    coordinates = [func if k == 0 else Derivative(func, (t,k))
                   for k in range(order, -1, -1) for func in funcs]
    rows, rhs = [], []
    for eq in eqs:
        result = linear_data(eq, coordinates, clear_denominators=False)
        if result is None or any(value.has(*funcs) for value in result[0]+[result[1]]):
            raise ODENonlinearError('The system of differential equations is nonlinear')
        coefficients, rest = result
        rows.append(coefficients)
        rhs.append(-rest)
    matrices = []
    for block in range(order+1):
        orientation = 1 if block == 0 else -1
        matrices.append(Matrix(len(eqs), len(funcs), lambda i,j:orientation*rows[i][block*len(funcs)+j]))
    return matrices, Matrix(len(eqs),1,rhs)


def canonical_odes(eqs, funcs, t):
    from ... import Derivative, Eq, Dummy, solve, cancel
    eqs, funcs, t = _inputs(eqs, funcs, t)
    if len(eqs) != len(funcs):
        raise NotImplementedError('Canonical systems require one equation per dependent function')
    orders = [max(ode_order(eq, func) for eq in eqs) for func in funcs]
    if max(orders, default=0) > 256:
        raise NotImplementedError('System order exceeds 256')
    derivatives = [Derivative(func,(t,order)) if order else func for func, order in zip(funcs, orders)]
    placeholders = [Dummy('highest_derivative') for _ in funcs]
    replaced = [eq.xreplace(dict(zip(derivatives,placeholders))) for eq in eqs]
    solutions = solve(replaced, placeholders, dict=True)
    if isinstance(solutions, dict):
        solutions = [solutions]
    output = []
    for solution in solutions:
        if not isinstance(solution, dict) or any(symbol not in solution for symbol in placeholders):
            raise NotImplementedError('The highest derivatives are not uniquely isolated in this branch')
        restore = dict(zip(placeholders,derivatives))
        output.append([Eq(derivative, solution[placeholder].xreplace(restore))
                       for derivative, placeholder in zip(derivatives,placeholders)])
    return output


def _zero_matrix(matrix):
    from ... import simplify, cancel
    return all(simplify(cancel(value)) == 0 for value in matrix)


def linodesolve(A, t, b=None, B=None, type='auto', doit=False, tau=None):
    from ... import Matrix, Dummy, integrate, Integral, sympify, S
    t = sympify(t)
    if not t.is_Symbol:
        raise ValueError('The independent coordinate must be a symbol')
    if tau is not None:
        raise NotImplementedError('Alternative system time coordinates are not implemented')
    if type not in ('auto','type1','type2','type3','type4'):
        raise NotImplementedError('This linear system solver type is not implemented')
    try:
        A = Matrix(A)
    except (ValueError, TypeError) as error:
        raise ValueError('A must be a concrete square matrix') from error
    if A.rows != A.cols:
        raise ValueError('A must be square')
    size = A.rows
    if size > 8:
        raise NotImplementedError('Matrix-exponential system solving is limited to eight functions')
    if b is None:
        b = Matrix(size,1,[S.Zero]*size)
    else:
        b = Matrix(b)
    if b.shape != (size,1):
        raise ValueError('The forcing must be a column vector matching A')
    constant = not any(value.has(t) for value in A)
    homogeneous = _zero_matrix(b)
    selected = ('type1' if homogeneous else 'type2') if constant else ('type3' if homogeneous else 'type4')
    if type != 'auto' and type != selected:
        raise ValueError('The requested solver type does not match the coefficient system')
    if B is None:
        B = A*t if constant else A.applyfunc(lambda value:integrate(value,t))
    else:
        B = Matrix(B)
        if B.shape != A.shape or not _zero_matrix(B.diff(t)-A):
            raise ValueError('B must be an antiderivative of A')
    if not _zero_matrix(A*B-B*A):
        raise NotImplementedError('A must commute with its matrix antiderivative')
    fundamental = B.exp()
    constants = Matrix([Dummy('C') for _ in range(size)])
    if homogeneous:
        result = fundamental*constants
    else:
        integrand = (-B).exp()*b
        primitive = integrand.applyfunc(lambda value:integrate(value,t) if doit else Integral(value,t))
        result = fundamental*(constants+primitive)
    return list(result)


def dsolve_system(eqs, funcs=None, ics=None):
    from ... import Eq, Symbol, Matrix, solve, default_sort_key, sympify
    eqs = [equation(eq) for eq in eqs]
    if not eqs:
        raise ValueError('Supply at least one differential equation')
    if funcs is None:
        funcs = sorted({node.expr for eq in eqs for node in nodes(eq) if node.is_Derivative
                        and node.expr.is_UndefinedFunction and len(node.expr.args)==1}, key=default_sort_key)
    funcs = list(funcs)
    if not funcs:
        raise ValueError('Cannot identify the dependent functions')
    t = funcs[0].args[0]
    eqs, funcs, t = _inputs(eqs, funcs, t)
    if len(eqs) != len(funcs):
        raise NotImplementedError('System solving requires a square differential system')
    (lead, coefficients), forcing = linear_ode_to_matrix(eqs, funcs, t, 1)
    inverse = lead.inv()
    rhs = linodesolve(inverse*coefficients, t, inverse*forcing, doit=True)
    occupied = {symbol.name for eq in eqs for symbol in eq.free_symbols}
    if ics:
        ics = {sympify(condition):sympify(value) for condition,value in ics.items()}
        occupied.update(symbol.name for pair in ics.items() for expression in pair
                        for symbol in expression.free_symbols)
    dummy_constants = sorted(set().union(*(value.free_symbols for value in rhs))
                             -set().union(*(eq.free_symbols for eq in eqs)), key=default_sort_key)
    replacements = {}
    index = 1
    for dummy in dummy_constants:
        while 'C'+str(index) in occupied:
            index += 1
        replacements[dummy] = Symbol('C'+str(index))
        occupied.add('C'+str(index))
        index += 1
    rhs = [value.xreplace(replacements) for value in rhs]
    if ics:
        constraints = []
        for condition, value in ics.items():
            index = next((i for i, func in enumerate(funcs)
                          if condition.func == func.func and len(condition.args)==1), None)
            if index is None:
                raise NotImplementedError('System initial conditions currently require function values')
            constraints.append(rhs[index].subs(t,condition.args[0])-value)
        values = solve(constraints,list(replacements.values()),dict=True)
        if isinstance(values,list) and len(values)==1:
            values = values[0]
        if not isinstance(values,dict):
            raise ValueError('Initial conditions do not select a unique supported solution branch')
        rhs = [value.subs(values) for value in rhs]
    return [Eq(func,value) for func,value in zip(funcs,rhs)]


__all__ = ['linear_ode_to_matrix','canonical_odes','linodesolve','ODEOrderError','ODENonlinearError']
