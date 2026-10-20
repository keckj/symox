"""Factorial-type combinatorial functions."""
from ... import (FallingFactorial, RisingFactorial, binomial, factorial, factorial2, ff, rf, subfactorial)

__all__ = ['FallingFactorial', 'RisingFactorial', 'binomial', 'factorial', 'factorial2', 'ff', 'rf', 'subfactorial']

from symox._function_heads import define_abstract as _define_abstract  # noqa: E402
_define_abstract(globals(), 'CombinatorialFunction MultiFactorial')
