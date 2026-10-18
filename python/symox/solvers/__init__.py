"""Native equation, system and solution-validation interfaces."""
from .solvers import *
from .solvers import __all__ as _scalar_names
from .solveset import solveset, solveset_real, solveset_complex, nonlinsolve
from .. import dsolve, checkodesol
from .decompogen import decompogen, compogen
from .deutils import ode_order, homogeneous_order
from .inequalities import (solve_poly_inequality,solve_rational_inequalities,
                           reduce_abs_inequality,reduce_abs_inequalities)

__all__=_scalar_names+['solveset','solveset_real','solveset_complex','nonlinsolve','dsolve','checkodesol','decompogen','compogen','ode_order','homogeneous_order',
                     'solve_poly_inequality','solve_rational_inequalities','reduce_abs_inequality','reduce_abs_inequalities']
