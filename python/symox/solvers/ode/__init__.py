"""Ordinary differential equations, following SymPy 1.14's sympy.solvers.ode.

The solver modules import the complete public surface, so the names resolve
on first access: the root import of symox loads this package before that
surface exists.
"""
from symox._api_exports import export as _export

_export(globals(), 'symox.solvers.ode.ode',
        'allhints checkinfsol classify_ode constantsimp dsolve homogeneous_order')
_export(globals(), 'symox.solvers.ode.lie_group', 'infinitesimals')
_export(globals(), 'symox.solvers.ode.subscheck', 'checkodesol')
_export(globals(), 'symox.solvers.ode.systems', 'canonical_odes linear_ode_to_matrix linodesolve')
_export(globals(), 'symox.solvers.ode.ode', 'classify_sysode solve_ics')
__all__ = [
    'allhints', 'checkinfsol', 'checkodesol', 'classify_ode', 'constantsimp',
    'dsolve', 'homogeneous_order', 'infinitesimals', 'canonical_odes', 'linear_ode_to_matrix',
    'linodesolve'
]
