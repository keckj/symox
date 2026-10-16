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
    result = _core._solve_poly_system(equations, variables)
    return None if result is None else [tuple(solution) for solution in result]


def nonlinsolve(system, *symbols):
    """Solve exact polynomial systems as a FiniteSet of solution tuples.

    Linear systems retain the existing parametric linsolve behavior. Nonlinear
    systems use the finite polynomial solver; positive-dimensional nonlinear
    varieties, transcendental equations, and incomplete radical root sets
    raise explicitly.
    """
    from . import EmptySet, FiniteSet, NonlinearError, _core
    from ._linear import linear_eq_to_matrix, linsolve

    variables = _symbols(symbols)
    if not variables:
        raise ValueError("Explicit symbols must be given")
    equations = _equations(system)
    if not equations:
        return EmptySet
    try:
        a, b = linear_eq_to_matrix(equations, *variables)
    except NonlinearError:
        result = _core._solve_poly_system(equations, variables)
        return EmptySet if result is None else FiniteSet(*(_core._apply("Tuple", *solution) for solution in result))
    return linsolve((a, b), *variables)


__all__ = ["solve_poly_system", "nonlinsolve"]
