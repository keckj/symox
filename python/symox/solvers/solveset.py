"""Set-valued solving through the native scalar equation engine."""
from .. import solveset, solveset_real, solveset_complex, linsolve, nonlinsolve

__all__=['solveset','solveset_real','solveset_complex','linsolve','nonlinsolve']

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox', 'NonlinearError linear_eq_to_matrix substitution')
