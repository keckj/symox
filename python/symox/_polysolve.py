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
    result = _core._nonlinsolve_parametric(equations, variables)
    rows = []
    for values, excluded in result:
        row = [_unwrap(Complement(Complexes, FiniteSet(*roots))) if roots else value
               for value, roots in zip(values, excluded)]
        rows.append(_core._apply("Tuple", *row))
    return EmptySet if not rows else FiniteSet(*rows)


__all__ = ["solve_poly_system", "nonlinsolve"]
