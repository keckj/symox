"""Scalar function heads applied directly to an abstract matrix expression.

SymPy 1.14 builds ``sin(A)`` for a ``MatrixSymbol`` ``A`` as an ordinary,
noncommutative function application, after the head's own ``eval`` rules
that apply to a matrix argument: odd heads take the sign out of ``-A``,
even heads drop it, ``erfc(-A)`` is ``2 - erfc(A)``, ``sign`` and ``arg``
drop a positive numeric factor. The application itself is the native one
(the matrix is held as a matrix proxy argument, as ``sin(x).subs(x, A)``
builds it). Heads whose reference evaluation rewrites a matrix argument
into another matrix expression (``Abs``, ``re``, ``im``) or reports a
reference error (``floor``, ``ceiling``, ``loggamma``, ``Heaviside``,
``DiracDelta``, ``Min``, ``Max``) are not admitted here.
"""

from . import _core

# Heads whose application to a matrix expression is built as in the reference.
ADMITTED = frozenset((
    'exp', 'log', 'sin', 'cos', 'tan', 'cot', 'sec', 'csc', 'asin', 'acos', 'atan', 'acot',
    'sinh', 'cosh', 'tanh', 'coth', 'asinh', 'acosh', 'atanh', 'gamma', 'erf', 'erfc',
    'sign', 'arg',
))
# f(-A) = -f(A) and f(-A) = f(A) in the reference evals.
_ODD = frozenset(('sin', 'tan', 'cot', 'csc', 'asin', 'atan', 'acot', 'sinh', 'tanh', 'coth',
                  'asinh', 'atanh', 'erf', 'sign'))
_EVEN = frozenset(('cos', 'sec', 'cosh'))
# Reference errors raised by the head's eval on a matrix argument.
_TYPE_ERRORS = frozenset(('floor', 'ceiling'))


def _minus_extractable(arg):
    """``could_extract_minus_sign`` of a matrix expression as in SymPy:
    a ``MatMul`` whose first argument is a negative number, and the
    reference rule of ``MatAdd``."""
    kind = type(arg).__name__
    if kind == 'MatMul':
        first = arg.args[0]
        return bool(getattr(first, 'is_Number', False) and first.is_extended_negative)
    if kind == 'MatAdd':
        return bool(arg.could_extract_minus_sign())
    return False


def _positive_factor(arg):
    """Replace the numeric leading factor ``c`` of a ``MatMul`` by its sign
    (``2*A`` -> ``A``, ``-2*A`` -> ``-A``), as ``sign.eval`` and ``arg.eval``
    do with ``as_coeff_Mul``. ``None`` when there is nothing to remove."""
    if type(arg).__name__ != 'MatMul':
        return None
    first = arg.args[0]
    if not getattr(first, 'is_Number', False) or first in (1, -1):
        return None
    if first.is_extended_positive or first.is_extended_negative:
        rest = arg.args[1:]
        from .matrices.expressions.matmul import MatMul
        inner = rest[0] if len(rest) == 1 else MatMul(*rest)
        return inner if first.is_extended_positive else -inner
    return None


def apply(name, arg, evaluate):
    """The application of head ``name`` to the matrix expression ``arg``,
    or raise the reference error. ``None`` when the head is not admitted."""
    if name in _TYPE_ERRORS:
        raise TypeError('Mix of Matrix and Scalar symbols')
    if name == 'loggamma':
        raise NotImplementedError()
    if name not in ADMITTED:
        return None
    if evaluate:
        if name in ('sec', 'csc') and type(arg).__name__ == 'ZeroMatrix':
            from .matrices.exceptions import NonInvertibleMatrixError
            raise NonInvertibleMatrixError('Matrix det == 0; not invertible.')
        if name in ('sign', 'arg'):
            inner = _positive_factor(arg)
            if inner is not None:
                return apply(name, inner, evaluate)
        if name in _ODD or name in _EVEN or name == 'erfc':
            if _minus_extractable(arg):
                value = apply(name, -arg, evaluate)
                if name in _ODD:
                    return -value
                if name == 'erfc':
                    return 2 - value
                return value
        return _core._apply(name, arg)
    return _core._apply(name, arg, evaluate=False)
