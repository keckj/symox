"""Trigonometric and inverse trigonometric functions
(sympy/functions/elementary/trigonometric.py)."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'sin cos tan cot sec csc sinc asin acos atan acot asec acsc atan2')

from symox._function_heads import define_abstract as _define_abstract  # noqa: E402
_define_abstract(globals(), 'TrigonometricFunction ReciprocalTrigonometricFunction '
                            'InverseTrigonometricFunction')
