"""Codec of the native ``convert_array_to_matrix`` passes and editor.

The passes run in Rust (``crates/symox-core/src/arrayexpr_a2m.rs``, binding
``_arrayexpr_a2m``). This module only encodes the expression tree, rebuilds
the objects from the result and, when the native port declines an input
(``_Decline`` while encoding, ``ArrayExprUnsupported`` from Rust), lets the
caller run the stock SymPy algorithm (``_stock``).
"""
import symox as S
from . import (ArrayAdd, ArrayContraction, ArrayDiagonal, ArrayElementwiseApplyFunc, ArraySymbol,
               ArrayTensorProduct, OneArray, PermuteDims, Reshape, ZeroArray, _shape)
from ._canon import _Codec
from .... import _core

_STRUCT = (ArrayTensorProduct, ArrayAdd, PermuteDims, ArrayContraction, ArrayDiagonal, ZeroArray)
Unsupported = _core.ArrayExprUnsupported


class _Decline(Exception):
    """An input the native port does not model."""


# A leaf the native port holds but no object here represents (SymPy's
# unprintable ``MatMul`` with an array coefficient, and its trace).
_HELD = object()


def _holds(n, leaves):
    if not isinstance(n, (tuple, list)):
        return False
    if len(n) == 2 and n[0] == 'L' and isinstance(n[1], int):
        return leaves[n[1]] is _HELD
    return any(_holds(x, leaves) for x in n)


def _as_lambda(f):
    if isinstance(f, _core.Expr) and getattr(f, 'func', None) is S.Lambda:
        return f
    d = S.Dummy('d')
    try:
        lam = S.Lambda(d, f(d))
    except Exception:
        raise _Decline
    if not isinstance(lam, _core.Expr):
        raise _Decline
    return lam


class _A2mCodec(_Codec):
    def __init__(self):
        super().__init__()
        self.descs = []
        self.opaque = {}

    def _leaf(self, desc, obj, key=None):
        if key is not None:
            k = self.opaque.get(key)
            if k is not None:
                return ('L', k)
            self.opaque[key] = len(self.descs)
        self.descs.append(desc)
        self.leaves.append(obj)
        return ('L', len(self.descs) - 1)

    def enc(self, e):
        t = type(e)
        if t in _STRUCT:
            return super().enc(e)
        if t is Reshape:
            return ('R', self.enc(e.expr), [self.dim(n) for n in e.shape])
        if t is ArrayElementwiseApplyFunc:
            f = e.function
            lam = _as_lambda(f)
            node = self.enc(e.expr)
            return ('F', node, self._leaf(('fn', lam, f), f)[1])
        if t is ArraySymbol:
            shape = _shape(e)
            return self._leaf(('a', e, len(shape), [self.dim(n) for n in shape]), e, ('a', e))
        if t is OneArray:
            return self._leaf(('one', [self.dim(n) for n in e.shape]), e)
        if isinstance(e, _core._MatrixExpr):
            return self._leaf(('m', e), e)
        if isinstance(e, S.Mul) and not e.is_commutative:
            return self._ncmul(e)
        if isinstance(e, bool):
            raise _Decline
        if isinstance(e, (int, float, complex)):
            e = S.sympify(e)
        if isinstance(e, _core.Expr):
            return self._leaf(('s', e), e)
        raise _Decline

    def _ncmul(self, e):
        # ``Mul(coeff, *arrays)``: the product SymPy forms of an array leaf
        # and a scalar (its array factors are the stock classes).
        from . import _stock
        coeff, ids = [], []
        for a in e.args:
            p = _stock.from_stock(a)
            if type(p) in (ArraySymbol, OneArray):
                ids.append(self.enc(p)[1])
            elif isinstance(a, _core.Expr) and a.is_commutative:
                coeff.append(a)
            else:
                raise _Decline
        self.descs.append(('ncmul', S.Mul(*coeff), ids))
        self.leaves.append(e)
        return ('L', len(self.descs) - 1)

    def _sync(self, descs, dims):
        from . import _stock
        self.dimlist = list(dims)
        for i in range(len(self.leaves), len(descs)):
            d = descs[i]
            if d[0] == 'one':
                self.leaves.append(OneArray(*[self.dimlist[k] for k in d[1]]))
            elif d[0] == 'held':
                self.leaves.append(_HELD)
            elif d[0] == 'ncmul':
                self.leaves.append(S.Mul(d[1], *[_stock.stock_leaf(self.leaves[k]) for k in d[2]]))
            else:
                self.leaves.append(d[1])

    def dec(self, n):
        if n is None:
            return None
        if n[0] == 'F':
            return ArrayElementwiseApplyFunc(self.leaves[n[2]], self.dec(n[1]))
        return super().dec(n)

    def order(self, nodes, descs=None, dims=None):
        if descs is not None:
            self._sync(descs, dims)
        if _holds(nodes, self.leaves):
            raise Unsupported('ordering an unevaluated matrix product')
        return super().order(nodes)

    def call(self, op, trees, extra=None):
        value, descs, dims = _core._arrayexpr_a2m(op, list(trees), self.descs, self.dimlist, self.order, extra)
        self._sync(descs, dims)
        return value


def run(op, exprs, extra=None, removed=False):
    """Runs pass ``op`` on ``exprs``; raises ``_Decline``/``Unsupported``."""
    c = _A2mCodec()
    trees = [c.enc(e) for e in exprs]
    value = c.call(op, trees, extra)
    if removed:
        return c.dec(value[0]), list(value[1])
    return c.dec(value)


def chain(expr):
    """``_chain_array_to_matrix`` (native ``A2m::chain``); raises
    ``_Decline``/``Unsupported`` outside the native port."""
    c = _A2mCodec()
    tree = c.enc(expr)
    symbols = []
    for k, leaf in enumerate(c.leaves):
        if type(leaf) is ArraySymbol and len(leaf.shape) == 2:
            try:
                symbols.append((k, S.MatrixSymbol(leaf.name, *leaf.shape)))
            except Exception:
                pass
    value = c.call('chain', [tree], symbols)
    if value == 'keep':
        return expr
    if value == 'fm':
        i, j = S.symbols('_i _j', integer=True)
        return S.FunctionMatrix(*expr.shape, S.Lambda((i, j), expr[i, j]))
    return c.dec(value)


def combine_removed(r1, r2):
    return list(_core._arrayexpr_a2m('combine_removed', [], [], [], None, (list(r1), list(r2))))
