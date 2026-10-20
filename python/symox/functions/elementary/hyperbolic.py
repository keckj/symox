"""Hyperbolic and inverse hyperbolic functions
(sympy/functions/elementary/hyperbolic.py)."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'sinh cosh tanh coth sech csch asinh acosh atanh acoth asech acsch')

from symox._function_heads import define_abstract as _define_abstract  # noqa: E402
_define_abstract(globals(), 'HyperbolicFunction ReciprocalHyperbolicFunction '
                            'InverseHyperbolicFunction')
