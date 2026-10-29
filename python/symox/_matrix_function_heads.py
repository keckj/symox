"""Scalar function heads applied directly to an abstract matrix expression.

SymPy 1.14 builds ``sin(A)`` for a ``MatrixSymbol`` ``A`` as an ordinary,
noncommutative function application, after the head's own ``eval`` rules
that apply to a matrix argument: odd heads take the sign out of ``-A``,
even heads drop it, ``erfc(-A)`` is ``2 - erfc(A)``, ``sign`` and ``arg``
drop a positive numeric factor. The application itself is the native one
(the matrix is held as a matrix proxy argument, as ``sin(x).subs(x, A)``
builds it). ``re``, ``im``, ``Abs``, ``Heaviside`` and ``DiracDelta``
go through ``conjugate(A)`` (``Adjoint(Transpose(A))``) as the reference
evals do; ``floor``, ``ceiling`` and ``loggamma`` report the reference
error.
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


# Heads rewritten through the conjugate of the matrix argument.
_CONJUGATE_HEADS = frozenset(('re', 'im', 'Abs', 'Heaviside', 'DiracDelta'))
# Matrix kinds equal to their own conjugate: Abs keeps them unevaluated.
_SELF_CONJUGATE = frozenset(('Identity', 'ZeroMatrix', 'OneMatrix'))
# Matrix kinds whose Abs is sqrt(expand_mul(M*conjugate(M))) with nothing
# to expand.
_ABS_LEAVES = frozenset(('MatrixSymbol', 'Transpose', 'Adjoint', 'MatPow', 'Inverse'))


def _is_leaf(arg):
    """A MatrixSymbol, possibly under Transpose, Adjoint, an integer MatPow
    or Inverse: the kinds whose ``M*conjugate(M)`` needs no expansion."""
    kind = type(arg).__name__
    if kind == 'MatrixSymbol':
        return True
    if kind in ('Transpose', 'Adjoint', 'Inverse'):
        return type(arg.args[0]).__name__ == 'MatrixSymbol'
    if kind == 'MatPow':
        return (type(arg.args[0]).__name__ == 'MatrixSymbol'
                and getattr(arg.args[1], 'is_Integer', False))
    return False


def _abs(arg):
    """``Abs.eval`` of a matrix expression as in SymPy 1.14, ``None`` when
    the argument is outside the admitted kinds."""
    kind = type(arg).__name__
    if kind in _SELF_CONJUGATE:
        # arg == conjugate(arg): Abs stays unevaluated.
        return _core._apply('Abs', arg, evaluate=False)
    if _is_leaf(arg):
        from .matrices.expressions.matmul import MatMul
        from . import S
        return MatMul(arg, arg.conjugate()).doit()**S.Half
    if kind != 'MatMul':
        return None
    from . import Abs, Mul, ImmutableMatrix
    known = []
    unknown = []
    for factor in arg.args:
        if getattr(factor, '_symox_matrixexpr', False):
            if type(factor).__name__ == 'MatPow' and not factor.args[1].is_positive:
                # The reference takes Pow(Abs(base), exp) here: not admitted.
                return None
            value = _abs(factor)
            if value is None or type(value).__name__ == 'Abs':
                return None
            known.append(value)
        elif getattr(factor, 'is_Matrix', False):
            known.append(ImmutableMatrix(factor.applyfunc(Abs)))
        else:
            # The reference MatMul keeps each scalar factor as an argument.
            for scalar in Mul.make_args(factor):
                value = Abs(scalar)
                if type(value).__name__ == 'Abs':
                    unknown.append(scalar)
                else:
                    known.append(value)
    result = known[0]
    for value in known[1:]:
        result = result*value
    if unknown:
        result = result*_core._apply('Abs', Mul(*unknown), evaluate=False)
    return result


def _real_check(name, arg):
    """Heaviside and DiracDelta of a matrix: the reference raises with the
    imaginary part of the argument."""
    imag = arg.as_real_imag()[1]
    message = ("Function defined only for Real Values. Complex part: %s  found in %s ."
               % (repr(imag), repr(arg)))
    if name == 'DiracDelta':
        from ._public_tools import filldedent
        message = filldedent(message)
    raise ValueError(message)


def conjugate_head(name, arg, evaluate):
    """re, im, Abs, Heaviside and DiracDelta of the matrix expression ``arg``
    (the reference rewrites them through ``conjugate(arg)``), ``None`` when
    the argument is not admitted."""
    if not evaluate:
        return _core._apply(name, arg, evaluate=False)
    if name in ('re', 'im'):
        return arg.as_real_imag()[0 if name == 're' else 1]
    if name == 'Abs':
        return _abs(arg)
    return _real_check(name, arg)


# re(f(M)) and im(f(M)) for a MatrixSymbol M (possibly transposed or
# adjoint), through f.as_real_imag in the reference.
_PART_MIX_ERRORS = frozenset(('sec', 'csc', 'cosh', 'tanh', 'coth', 'acosh', 'atanh'))


def function_part(part, value):
    """``re``/``im`` of a scalar function of a matrix, or ``None``."""
    head = getattr(getattr(value, 'func', None), '__name__', None)
    args = getattr(value, 'args', ())
    if head is None or len(args) != 1:
        return None
    inner = args[0]
    kind = type(inner).__name__
    if not (kind == 'MatrixSymbol' or kind in ('Transpose', 'Adjoint')
            and type(inner.args[0]).__name__ == 'MatrixSymbol'):
        return None
    if head == 'exp':
        raise TypeError('MatAdd.expand() takes 1 positional argument but 2 were given')
    if head in _PART_MIX_ERRORS:
        raise TypeError('Mix of Matrix and Scalar symbols')
    from . import S
    if head == 'log':
        if part == 're':
            return _core._apply('log', _abs(inner))
        return _core._apply('arg', inner)
    if head == 'arg':
        return value if part == 're' else S.Zero
    return None


def apply(name, arg, evaluate):
    """The application of head ``name`` to the matrix expression ``arg``,
    or raise the reference error. ``None`` when the head is not admitted."""
    if name in _CONJUGATE_HEADS:
        return conjugate_head(name, arg, evaluate)
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
