"""Scalar solvers and solution checks using native kernels."""
from .. import solve, nsolve
from .._solver_checks import checksol, solve_linear, solve_linear_system_LU, solve_undetermined_coeffs
from .._linear import solve_linear_system, linear_eq_to_matrix, linsolve

__all__=['solve','nsolve','checksol','solve_linear','solve_linear_system',
    'solve_linear_system_LU','linear_eq_to_matrix','linsolve','solve_undetermined_coeffs']
