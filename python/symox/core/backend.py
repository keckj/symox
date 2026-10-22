"""SymPy's ``sympy.core.backend``: the core names used by ``physics.mechanics``.

The reference switches to SymEngine under ``USE_SYMENGINE``; symox always
provides its own objects.
"""
from symox import (Add, Basic, diff, Function, expand, Derivative, Mul, igcd,
                   pi, I, Integer, Rational, E, S, Symbol, var, symbols,
                   SympifyError, sympify, log, exp, coth, sinh, acosh, acoth,
                   tanh, asinh, atanh, cosh, sqrt, csc, asec, cos, atan, sec,
                   acot, asin, tan, sin, cot, acsc, acos, gamma, eye, zeros,
                   diag, Matrix, ones, symarray, ImmutableMatrix, lambdify)
from symox.core.function import AppliedUndef
from symox.matrices import MatrixBase

USE_SYMENGINE = False


def _simplify_matrix(M):
    """Return a simplified copy of the matrix M"""
    if not isinstance(M, (Matrix, ImmutableMatrix)):
        raise TypeError("The matrix M must be an instance of Matrix or ImmutableMatrix")
    Mnew = M.as_mutable()
    Mnew.simplify()
    if isinstance(M, ImmutableMatrix):
        Mnew = Mnew.as_immutable()
    return Mnew
