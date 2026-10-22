"""Polynomial system interfaces backed by native elimination and radical roots."""

from ._linear import _equations, _symbols


def solve_poly_system(seq, *gens, strict=False, **options):
    """Return the complete finite list of polynomial solution tuples.

    Exact rational systems use native lexicographic Groebner elimination.
    Symbolic coefficients support linear elimination and triangular root
    systems. Inconsistent systems return None. Positive-dimensional systems
    and roots outside the implemented radical algorithms raise explicitly;
    partial root lists are never returned, including when strict is False.
    """
    from . import _core, default_sort_key

    if options:
        raise NotImplementedError("Polynomial system domain and extension options are not implemented")
    if not isinstance(strict, bool):
        raise TypeError("strict must be a bool")
    equations = _equations(seq)
    variables = _symbols(gens)
    if not variables:
        variables = tuple(sorted({v for e in equations for v in e.free_symbols}, key=default_sort_key))
    if not variables:
        raise ValueError("Polynomial system solving requires generators")
    try:
        result = _core._solve_poly_system(equations, variables)
    except NotImplementedError as error:
        if strict and 'cannot be expressed by the implemented radical root' in str(error):
            # SymPy's strict mode reports roots(strict=True)'s error class.
            from .polys.polyerrors import UnsolvableFactorError
            raise UnsolvableFactorError(
                '\nStrict mode: some factors cannot be solved in radicals, so a complete\n'
                'list of solutions cannot be returned. Call roots with strict=False to\n'
                'get solutions expressible in radicals (if there are any).') from error
        raise
    return None if result is None else [tuple(solution) for solution in result]


def nonlinsolve(system, *symbols):
    """Solve supported exact rational systems, including free parameters.

    Native factorization and triangular linear/quadratic elimination retain
    zero-coefficient branches. Free parameters retain their symbols; excluded
    parameter values are represented by Complement(Complexes, FiniteSet(...)).
    Every branch must be solved completely, otherwise NotImplementedError is
    raised. Conditions on external parameters require including those symbols
    among the unknowns.
    """
    from . import Complexes, Complement, EmptySet, FiniteSet, _core, _unwrap

    variables = _symbols(symbols)
    if not variables:
        raise ValueError("Explicit symbols must be given")
    equations = _equations(system)
    if not equations:
        return EmptySet
    if any(not getattr(v, 'is_Symbol', False) for v in variables):
        # SymPy recasts the unknowns that are not Symbols (f(x), x + 2, ...)
        # to Dummy symbols and restores them in the solution tuples.
        from .solvers._algebraic import recast_to_symbols
        recast, unknowns, swap = recast_to_symbols(equations, variables)
        solutions = nonlinsolve(recast, unknowns)
        if not isinstance(solutions, FiniteSet):
            return solutions.xreplace(swap)
        return FiniteSet(*[_core._apply("Tuple", *(i.xreplace(swap) for i in s)) for s in solutions])
    if any(_core._has_noncommutative_atoms(e) for e in equations):
        return _noncommutative_nonlinsolve(equations, variables)
    try:
        result = _core._nonlinsolve_parametric(equations, variables)
    except (ValueError, NotImplementedError):
        # A system with a non-polynomial equation (exp, log, trigonometric
        # functions of the unknowns) takes the reference branch: Groebner
        # elimination of the polynomial part, then substitution, which keeps
        # the periodic ImageSet families.
        from .solvers._substitution import _is_polynomial_system, nonlinsolve_substitution
        if _is_polynomial_system(equations, variables):
            solutions = _zero_dimensional_numerators(equations, variables)
            if solutions is None:
                raise
            return solutions
        return nonlinsolve_substitution(list(equations), list(variables))
    rows = []
    for values, excluded in result:
        row = [_unwrap(Complement(Complexes, FiniteSet(*roots))) if roots else value
               for value, roots in zip(values, excluded)]
        rows.append(_core._apply("Tuple", *row))
    return EmptySet if not rows else FiniteSet(*rows)


def _zero_dimensional_numerators(equations, variables):
    """The reference zero-dimensional branch of ``nonlinsolve`` for a
    system whose numerators are polynomials with exact numeric coefficients:
    grevlex Groebner basis of the numerators, FGLM to lex, then the complete
    ``solve_poly_system`` roots. Unlike the reference, roots where a
    denominator of the original system vanishes are dropped (they are not
    solutions). Returns None outside that class, so the caller raises."""
    from . import EmptySet, FiniteSet, _core, groebner, simplify, solve_poly_system, together

    numerators, denominators = [], []
    for eq in equations:
        if getattr(eq, 'is_Equality', False):
            eq = eq.lhs - eq.rhs
        numerator, denominator = together(eq).as_numer_denom()
        if not numerator.free_symbols <= set(variables) or not denominator.free_symbols <= set(variables):
            return None
        numerators.append(numerator)
        if denominator.free_symbols:
            denominators.append(denominator)
    try:
        basis = groebner(numerators, *variables, order='grevlex')
        if 1 in list(basis):
            return EmptySet
        if not basis.is_zero_dimensional:
            return None
        roots = solve_poly_system(list(basis.fglm('lex')), *variables)
    except (NotImplementedError, ValueError):
        return None
    if roots is None:
        return EmptySet
    rows = []
    for root in roots:
        values = dict(zip(variables, root))
        vanishing = [simplify(d.subs(values)).is_zero for d in denominators]
        if any(v is None for v in vanishing):
            return None
        if not any(vanishing):
            rows.append(_core._apply("Tuple", *root))
    return FiniteSet(*rows) if rows else EmptySet


def _noncommutative_nonlinsolve(equations, variables):
    """One equation with noncommutative coefficients in one unknown, solved
    by the noncommutative ``solve`` (commutative unknowns take commutative
    values, a noncommutative unknown keeps the factor order). An empty
    result is a proven absence of solutions; other systems raise."""
    from . import EmptySet, FiniteSet, _core, solve

    if len(equations) != 1 or len(variables) != 1:
        raise NotImplementedError(
            "systems with noncommutative coefficients are not implemented")
    solutions = solve(equations[0], variables[0], dict=False)
    if not solutions:
        return EmptySet
    return FiniteSet(*[_core._apply("Tuple", s) for s in solutions])


__all__ = ["solve_poly_system", "nonlinsolve"]
