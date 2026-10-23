"""SymPy linear-system interfaces backed by native coefficient extraction/RREF."""

import itertools


def _symbols(values, count=None):
    from . import DenseMatrix as Matrix, sympify
    from ._core import Expr

    if len(values) == 1 and not isinstance(values[0], Expr):
        source = values[0]
        if isinstance(source, (set, frozenset, dict)):
            raise TypeError("Symbols must be an ordered sequence")
        if isinstance(source, (list, tuple, Matrix)):
            values = tuple(source)
        elif count is not None and hasattr(source, "__next__"):
            values = tuple(itertools.islice(source, count))
    result = tuple(sympify(value) for value in values)
    if len(set(result)) != len(result):
        raise ValueError("Symbols must be unique")
    return result


def _equations(value):
    from . import DenseMatrix as Matrix, sympify
    from ._core import Expr

    if isinstance(value, (Expr, str, int, float)):
        return [sympify(value)]
    if isinstance(value, (set, frozenset, dict)):
        raise TypeError("Equations must be an ordered sequence")
    return [sympify(item) for item in value]


def linear_eq_to_matrix(equations, *symbols):
    """Extract ``A, b`` from equations ``A*x=b`` without nonlinear cancellation.

    Coefficients can be symbolic expressions independent of the requested
    unknowns. Nonlinear terms raise ``NonlinearError``.
    """
    from . import DenseMatrix as Matrix, _core

    variables = _symbols(symbols)
    if not variables:
        raise ValueError("Symbols must be given")
    expressions = _equations(equations)
    try:
        coefficients, constants = _core._linear_eq_coefficients(expressions, variables)
    except _core.NonlinearError as error:
        message = _nonlinear_message(expressions, set(variables))
        if message is None:
            raise
        raise _core.NonlinearError(message) from None
    return Matrix(len(expressions), len(variables), coefficients), Matrix(len(expressions), 1, constants)


class _Nonlinear(Exception):
    pass


def _lin_eq2dict(a, symset):
    # SymPy's polys.matrices.linsolve._lin_eq2dict, used only to report the
    # same offending term as the reference once nonlinearity is detected.
    if a in symset:
        return True
    if a.is_Add:
        return any([_lin_eq2dict(ai, symset) for ai in a.args])
    if a.is_Mul:
        found = False
        for ai in a.args:
            if _lin_eq2dict(ai, symset):
                if found:
                    raise _Nonlinear('\nnonlinear cross-term: %s' % a)
                found = True
        return found
    if not a.has_xfree(symset):
        return False
    raise _Nonlinear('nonlinear term: %s' % a)


def _nonlinear_message(expressions, symset):
    try:
        for e in expressions:
            for side in ((e.lhs, e.rhs) if getattr(e, 'is_Equality', False) else (e,)):
                _lin_eq2dict(side, symset)
    except _Nonlinear as error:
        return str(error)
    except Exception:
        return None
    return None


def _matrix_system(system):
    from . import DenseMatrix as Matrix

    if isinstance(system, Matrix) or (hasattr(system, "rows") and hasattr(system, "cols")):
        augmented = Matrix(system)
        if augmented.cols == 0:
            raise ValueError("An augmented matrix must include a right-hand-side column")
        return augmented[:, :-1], augmented[:, -1:]
    if isinstance(system, (tuple, list)) and len(system) == 2:
        a, b = system
        if all(isinstance(item, Matrix) or (hasattr(item, "rows") and hasattr(item, "cols")) for item in (a, b)):
            return Matrix(a), Matrix(b)
    return None


def _parameters(a, b, count):
    from . import Symbol

    names = {symbol.name for expression in itertools.chain(a, b) for symbol in expression.free_symbols}
    base = "tau"
    while any(name.startswith(base) and name[len(base):].isdigit() for name in names):
        base = "_" + base
    return tuple(Symbol(f"{base}{index}") for index in range(count))


def linsolve(system, *symbols):
    """Solve exact linear systems, including inconsistent and parametric ones.

    Accepts an equation sequence with explicit unknowns, an augmented Matrix,
    or ``(A, b)`` matrices. Returns a FiniteSet of a solution tuple or EmptySet.
    Free unknowns retain the supplied symbols; matrix-only calls generate tau
    parameters without colliding with coefficient symbols.
    """
    from . import EmptySet, FiniteSet, _core

    matrices = _matrix_system(system)
    if matrices is None:
        expressions = _equations(system)
        if not expressions:
            return EmptySet
        variables = _symbols(symbols)
        if not variables:
            raise ValueError("Explicit symbols must be given when passing equations")
        from . import _noncommutative
        if _noncommutative._any_nc(expressions):
            # linear_eq_to_matrix -> Poly: the reference rejects nc input.
            raise _noncommutative._polynomial_error()
        a, b = linear_eq_to_matrix(expressions, *variables)
    else:
        a, b = matrices
        if b.cols != 1 or a.rows != b.rows:
            raise ValueError("The right-hand side must be a matching column matrix")
        variables = _symbols(symbols, a.cols)
        if variables and len(variables) < a.cols:
            raise ValueError("At least as many symbols as columns of A must be given")
        if variables:
            variables = variables[:a.cols]
        else:
            # The native RREF uses one parameter per free column. First solve
            # with distinct placeholders, then compact only the free names.
            variables = _parameters(a, b, a.cols)
            solution = _core._linsolve_matrix(a, b, variables)
            if solution is None:
                return EmptySet
            free = [variable for variable, value in zip(variables, solution) if value == variable]
            parameters = _parameters(a, b, len(free))
            mapping = dict(zip(free, parameters))
            solution = [value.xreplace(mapping) for value in solution]
            return FiniteSet(_core._apply("Tuple", *solution))
    solution = _core._linsolve_matrix(a, b, variables)
    return EmptySet if solution is None else FiniteSet(_core._apply("Tuple", *solution))


def solve_linear_system(system, *symbols, **flags):
    """Solve an augmented linear Matrix; free unknowns are omitted from the dict."""
    from . import DenseMatrix as Matrix, _core

    if flags:
        raise NotImplementedError("solve_linear_system options are not implemented")
    matrix = Matrix(system)
    variables = _symbols(symbols)
    if matrix.cols != len(variables)+1:
        raise ValueError("Augmented matrix columns must match the symbols plus one")
    solution = _core._linsolve_matrix(matrix[:, :-1], matrix[:, -1:], variables)
    if solution is None:
        return None
    return {variable: value for variable, value in zip(variables, solution) if value != variable}


__all__ = ["linear_eq_to_matrix", "linsolve", "solve_linear_system"]
