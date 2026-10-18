"""Native equation, system and solution-validation interfaces."""
from .solvers import *
from .solvers import __all__ as _scalar_names
from .solveset import solveset, solveset_real, solveset_complex, nonlinsolve
from .. import dsolve, checkodesol
from .decompogen import decompogen, compogen
from .deutils import ode_order, homogeneous_order
from .ode import classify_ode, dsolve, checkodesol
from .pde import pdsolve, classify_pde, checkpdesol, pde_separate, pde_separate_add, pde_separate_mul
from .inequalities import (solve_poly_inequality,solve_rational_inequalities,
                           reduce_abs_inequality,reduce_abs_inequalities)

__all__=_scalar_names+['solveset','solveset_real','solveset_complex','nonlinsolve','dsolve','checkodesol','decompogen','compogen','ode_order','homogeneous_order',
                     'solve_poly_inequality','solve_rational_inequalities','reduce_abs_inequality','reduce_abs_inequalities']
__all__ += ['classify_ode', 'pdsolve', 'classify_pde', 'checkpdesol', 'pde_separate',
            'pde_separate_add', 'pde_separate_mul']

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox', 'Complexes check_assumptions diophantine factor_system failing_assumptions reduce_inequalities rsolve rsolve_hyper rsolve_poly rsolve_ratio solve_poly_system solve_triangulated solve_univariate_inequality substitution')

# Load the homonymous package before binding the solver function.
from .diophantine import diophantine
from .simplex import linprog, lpmin, lpmax
__all__ += ['linprog', 'lpmin', 'lpmax']
