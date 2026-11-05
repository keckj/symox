"""Stock fallback of the array-expression conversions.

``symox._arrayexpr`` is the unmodified SymPy 1.14 implementation, built on
its own (``Basic``-based) array-expression classes. When the native (Rust)
conversion declines an input (an object kind it does not model), the public
expression is translated structurally into those classes, the stock
algorithm runs, and the result is translated back. Nothing is canonicalized
on the way: the translation keeps every node as given.
"""
from . import (ArrayAdd, ArrayContraction, ArrayDiagonal, ArrayElement, ArrayElementwiseApplyFunc, ArraySymbol,
               ArrayTensorProduct, OneArray, PermuteDims, Reshape, ZeroArray, _ArrayExpr)


def _v():
    from ...._arrayexpr import array_expressions as V
    return V


def stock_leaf(e):
    """``symox._arrayexpr`` class of an ``ArraySymbol`` or ``OneArray``: the
    factors of the ``Mul`` SymPy forms of an array leaf and a scalar."""
    V = _v()
    if type(e) is ArraySymbol:
        return V.ArraySymbol(e.name, e.shape)
    return V.OneArray(*e.shape)


def to_stock(e):
    """Public array expression -> ``symox._arrayexpr`` expression."""
    V = _v()
    t = type(e)
    if t is ArraySymbol or t is OneArray:
        return stock_leaf(e)
    if t is ZeroArray:
        return V.ZeroArray(*e.shape)
    if t is ArrayTensorProduct:
        return V.ArrayTensorProduct(*[to_stock(a) for a in e.args])
    if t is ArrayAdd:
        return V.ArrayAdd(*[to_stock(a) for a in e.args])
    if t is PermuteDims:
        return V.PermuteDims(to_stock(e.expr), e.permutation)
    if t is ArrayContraction:
        return V.ArrayContraction(to_stock(e.expr), *e.groups)
    if t is ArrayDiagonal:
        return V.ArrayDiagonal(to_stock(e.expr), *e.groups)
    if t is Reshape:
        return V.Reshape(to_stock(e.expr), e.shape)
    if t is ArrayElementwiseApplyFunc:
        return V.ArrayElementwiseApplyFunc(e.function, to_stock(e.expr))
    if t is ArrayElement:
        return V.ArrayElement(to_stock(e.name), e.indices)
    if isinstance(e, (list, tuple)):
        return type(e)(to_stock(a) for a in e)
    return e


def from_stock(e):
    """``symox._arrayexpr`` expression -> public array expression."""
    V = _v()
    t = type(e)
    if t is V.ArraySymbol:
        return ArraySymbol(e.name, e.shape)
    if t is V.ZeroArray:
        return ZeroArray(*e.shape)
    if t is V.OneArray:
        return OneArray(*e.shape)
    if t is V.ArrayTensorProduct:
        return ArrayTensorProduct(*[from_stock(a) for a in e.args])
    if t is V.ArrayAdd:
        return ArrayAdd(*[from_stock(a) for a in e.args])
    if t is V.PermuteDims:
        return PermuteDims(from_stock(e.expr), e.permutation)
    if t is V.ArrayContraction:
        return ArrayContraction(from_stock(e.expr), *e.contraction_indices)
    if t is V.ArrayDiagonal:
        return ArrayDiagonal(from_stock(e.expr), *e.diagonal_indices)
    if t is V.Reshape:
        return Reshape(from_stock(e.expr), e.shape)
    if t is V.ArrayElementwiseApplyFunc:
        return ArrayElementwiseApplyFunc(e.function, from_stock(e.expr))
    if t is V.ArrayElement:
        return ArrayElement(from_stock(e.name), e.indices)
    if isinstance(e, (list, tuple)):
        return type(e)(from_stock(a) for a in e)
    return e


def stock_call(module, name, *args):
    """Runs ``symox._arrayexpr.<module>.<name>`` on translated arguments."""
    import importlib
    mod = importlib.import_module('symox._arrayexpr.' + module)
    return from_stock(getattr(mod, name)(*[to_stock(a) for a in args]))


def is_public(e):
    return isinstance(e, _ArrayExpr) or isinstance(e, ArrayElement)
