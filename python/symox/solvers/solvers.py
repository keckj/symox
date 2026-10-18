"""Scalar solvers and solution checks using native kernels."""
from .. import solve, nsolve
from .._solver_checks import checksol, solve_linear, solve_linear_system_LU, solve_undetermined_coeffs
from .._linear import solve_linear_system, linear_eq_to_matrix, linsolve
from ._algebraic import (denoms,recast_to_symbols,det_perm,det_minor,det_quick,
                         inv_quick,minsolve_linear_system,unrad)

__all__=['solve','nsolve','checksol','solve_linear','solve_linear_system',
    'solve_linear_system_LU','linear_eq_to_matrix','linsolve','solve_undetermined_coeffs']
__all__ += ['denoms','recast_to_symbols','det_perm','det_minor','det_quick','inv_quick',
            'minsolve_linear_system','unrad']
