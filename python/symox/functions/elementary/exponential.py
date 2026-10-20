"""Exponential, logarithm and Lambert W functions."""
from ... import (LambertW, exp, exp_polar, log)

__all__ = ['LambertW', 'exp', 'exp_polar', 'log']

from symox._function_heads import define_abstract as _define_abstract  # noqa: E402
_define_abstract(globals(), 'ExpBase')
