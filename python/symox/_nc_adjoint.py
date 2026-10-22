"""Adjoint and transpose of noncommutative scalar expressions.

A port of SymPy's ``adjoint`` and ``transpose`` classes together with the
``_eval_adjoint``/``_eval_transpose`` rules of ``Expr``, ``Add``, ``Mul``,
``Pow``, ``exp``, ``Piecewise`` and ``conjugate``: products are reversed,
sums are mapped termwise, and an atom without a rule stays as an
unevaluated ``adjoint(A)``/``transpose(A)`` node. Commutative expressions
never reach this module; their adjoint is the native conjugate and their
transpose the expression itself.
"""
import symox as _sp
from symox import _core


def _has_nc(e):
    return isinstance(e, _core.Expr) and _core._has_noncommutative_atoms(e)


def _native_conjugate(e):
    return _core._apply('conjugate', e)


def _is_head(e, cls):
    return type(e) is cls or getattr(e, 'func', None) is cls


def _conjugate_arg(e):
    """The argument of an unevaluated conjugate node, else None."""
    if getattr(e, 'is_Function', False) and type(e).__name__ != 'adjoint' \
            and type(e).__name__ != 'transpose':
        try:
            name = e.func.__name__
        except AttributeError:
            return None
        if name == 'conjugate' and len(e.args) == 1:
            return e.args[0]
    return None


def adjoint_of(e):
    """SymPy's ``adjoint(e)`` for a noncommutative ``e``."""
    if not _has_nc(e):
        return _native_conjugate(e)
    return adjoint(e)


def transpose_of(e):
    """SymPy's ``transpose(e)`` for a noncommutative ``e``."""
    if not _has_nc(e):
        return e
    return transpose(e)


def _eval_adjoint(e):
    # Expr._eval_adjoint and the overrides of the node classes
    if e.is_Mul:
        return _sp.Mul(*[adjoint_of(t) for t in reversed(e.args)])
    if e.is_Add:
        return _sp.Add(*[adjoint_of(t) for t in e.args])
    if e.is_Pow:
        base, exponent = e.args
        if exponent.is_integer:
            return adjoint_of(base)**exponent
        if base.is_positive:
            return base**adjoint_of(exponent)
        return None
    if isinstance(e, adjoint):
        return e.args[0]
    if isinstance(e, transpose):
        return _native_conjugate(e.args[0])
    arg = _conjugate_arg(e)
    if arg is not None:
        return transpose_of(arg)
    if getattr(e, 'is_Piecewise', False):
        return _sp.Piecewise(*[(adjoint_of(x), c) for x, c in e.args])
    if getattr(e, 'func', None) is _sp.exp:
        return _sp.exp(adjoint_of(e.args[0]))
    if e.is_hermitian:
        return e
    if e.is_antihermitian:
        return -e
    # Expr._eval_conjugate (real or imaginary) cannot hold for a
    # noncommutative expression; Expr._eval_transpose then needs
    # hermiticity, already tested above.
    return None


def _eval_transpose(e):
    if e.is_Mul:
        return _sp.Mul(*[transpose_of(t) for t in reversed(e.args)])
    if e.is_Add:
        return _sp.Add(*[transpose_of(t) for t in e.args])
    if e.is_Pow:
        base, exponent = e.args
        if base.is_complex or base.is_infinite:
            return e
        if exponent.is_integer:
            return transpose_of(base)**exponent
        return None
    if isinstance(e, adjoint):
        return _native_conjugate(e.args[0])
    if isinstance(e, transpose):
        return e.args[0]
    arg = _conjugate_arg(e)
    if arg is not None:
        return adjoint_of(arg)
    if getattr(e, 'is_Piecewise', False):
        return _sp.Piecewise(*[(transpose_of(x), c) for x, c in e.args])
    if getattr(e, 'func', None) is _sp.exp:
        return _sp.exp(transpose_of(e.args[0]))
    if e.is_hermitian:
        return conjugate_of(e)
    if e.is_antihermitian:
        return -conjugate_of(e)
    return None


def conjugate_of(e):
    """``conjugate(e)`` with the conjugate rules of the adjoint and
    transpose nodes (``conjugate(adjoint(A)) == transpose(A)``)."""
    value = _native_conjugate(e)
    if not _has_nc(value):
        return value
    return _fold_conjugates(value)


def _fold_conjugates(e):
    if not isinstance(e, _core.Expr) or not e.args:
        return e
    arg = _conjugate_arg(e)
    if arg is not None:
        if isinstance(arg, adjoint):
            return transpose_of(arg.args[0])
        if isinstance(arg, transpose):
            return adjoint_of(arg.args[0])
        return e
    args = [_fold_conjugates(a) for a in e.args]
    if all(a is b or a == b for a, b in zip(args, e.args)):
        return e
    return e.func(*args)


class adjoint(_sp.Function):
    """Conjugate transpose of a noncommutative scalar (unevaluated node)."""
    is_commutative = False
    _symox_user_hooks = True

    @classmethod
    def eval(cls, arg):
        if not _has_nc(arg):
            return _native_conjugate(arg)
        obj = _eval_adjoint(arg)
        if obj is not None:
            return obj
        obj = _eval_transpose(arg)
        if obj is not None:
            return _native_conjugate(obj)

    def _eval_subs(self, old, new):
        if self == old:
            return new
        return adjoint_of(self.args[0].subs(old, new))


class transpose(_sp.Function):
    """Transpose of a noncommutative scalar (unevaluated node)."""
    is_commutative = False
    _symox_user_hooks = True

    @classmethod
    def eval(cls, arg):
        if not _has_nc(arg):
            return arg
        obj = _eval_transpose(arg)
        if obj is not None:
            return obj

    def _eval_subs(self, old, new):
        if self == old:
            return new
        return transpose_of(self.args[0].subs(old, new))

